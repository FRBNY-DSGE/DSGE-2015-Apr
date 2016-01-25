% OVERVIEW
%
% This is a dsge likelihood function that can handle 2-part estimation where
% there is a model switch.
%
% HOWEVER, this program is by and large the same as objfcnmhdsge_2part.m, save
% the bound-checking. Therefore, the code has been substantially consolidated
% and much of the code supporting the operations in this function are shared by
% objfcnmhdsge_2part.m
%
% Note: given the multi-period setup, all relevant matrices will be stored in
% an array structure since they will change over time. As currently
% implemented, there are 3 time periods, which index the structure: presample,
% normal sample, ZB sample.
%
% If there is no model switch, then we filter over the main sample all at once.
%
function pyt = dsgelh(para,YY,YY0,nobs,nlags,nvar,mspec,npara,coint,cointadd,YYcoint0,nant,antlags)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: Set up structures that will hold all of the transition equation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the general classes of 2-part estimation models where
% operations/indexing is similar within a particular class.
class2part;

% Create structure to hold all matrices -- data matrices, state equation matrices, measurement matrices, etc.
[mt, pd] = dsgelh_partition(YY0, YY, nvar, nant, antlags);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
%% transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the model
[mt(pd).TTT, mt(pd).RRR, mt(pd).CCC, valid] = dsgesolv(mspec,para, mt(end).nant{:});

%% For models using 2part estimation: Get the normal, no ZB model matrices
if any(mspec == class2part_all)

  % Get the starting indices
  [start_ant_state, start_ant_shock, revol_ind] = get_start_ant(mspec, nant);

  [mt(pd-1).TTT, mt(pd-1).RRR, mt(pd-1).CCC]  = ...
    dsgelh_getNoZB(nant, start_ant_state, start_ant_shock, mt(pd).TTT, mt(pd).RRR, mt(pd).CCC,revol_ind{:});
end

if valid < 1;
    pyt = -1E10;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some helper functions to create compound matrices by passing the
% period index to access the relevant EE, MM, QQ, RR, VV matrices
makeHH = @(E, M, Q) E + M*Q*M';
makeVV = @(Q, M) Q*M';
makeVVall = @(R, Q, V, H) [[R*Q*R', R*V]; ...
                           [V'*R', H]];


% Get measurement equation matrices set up for all periods that aren't the presample
for p = 2:pd

  % Save measurement equations
  [mt(p).ZZ, mt(p).DD, mt(p).DDcointadd, mt(p).QQ, mt(p).EE, mt(p).MM, retcode] = ...
    feval(['measur',num2str(mspec)], mt(p).TTT, mt(p).RRR,valid,para,mt(p).nvar,nlags,mspec,npara,coint,cointadd,mt(p).nant{:});

  if retcode == 0
      % invalid parameterization
      pyt = -1E10;
      return;
  end;

  mt(p).HH = makeHH(mt(p).EE, mt(p).MM, mt(p).QQ);
  mt(p).VV = makeVV(mt(p).QQ, mt(p).MM);
  mt(p).VVall = makeVVall(mt(p).RRR, mt(p).QQ, mt(p).VV, mt(p).HH);

  if p == 2 && any(mt(p).CCC ~= 0)
    mt(p).DD = mt(p).DD + (mt(p).ZZ)*((eye(size(mt(p).TTT))-mt(p).TTT)\mt(p).CCC);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: compute log-likelihood using Kalman filter - written by Iskander
%%         note that Iskander's program assumes a transition equation written as:
%%         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
%%         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
%%         and  VV2 = cov(eps2_t,u_u) = RRR*VV
%%         define VVall as the joint variance of the two shocks VVall = var([eps2_t;u_t])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PRESAMPLE

% Solve lyapunov with normal period state matrices (i.e. period 2 matrices)
[A0,P0] = lyap_nonstationary(mspec,para,mt(2).TTT,mt(2).RRR,mt(2).QQ);

% If there is a presample
if ~isempty(mt(1).YY)

  % If number of presample series differ from number of measurement equation
  % series, chop off extra measurement equations and recompute HH, VV, VVall
  % for presample. Otherwise, HH, VV, VVall same in presample as in normal pd
  if size(mt(1).YY,2) ~= size(mt(2).ZZ,1)
    nvar0 = size(mt(1).YY,2);

    % Chop off parts from measurement matrices
    meas_mats = {'ZZ', 'DD', 'EE', 'MM'};
    for mm = 1:length(meas_mats)
      mt(1).(meas_mats{mm}) = mt(2).(meas_mats{mm})(1:nvar0,:);
    end
    mt(1).EE = mt(1).EE(:,1:nvar0);

    % Recompute
    mt(1).HH = makeHH(mt(1).EE, mt(1).MM, mt(2).QQ);
    mt(1).VV = makeVV(mt(2).QQ, mt(1).MM);
    mt(1).VVall = makeVVall(mt(2).RRR, mt(2).QQ, mt(1).VV, mt(1).HH);

  else
    mt(1).ZZ = mt(2).ZZ;
    mt(1).DD = mt(2).DD;
    mt(1).VVall = mt(2).VVall;
  end

  %% run Kalman filter on initial observations
  [mt(1).pyt0, mt(1).zend, mt(1).Pend] = ...
        kalcvf2NaN( (mt(1).YY)', 1, zeros(size(mt(2).TTT,2),1), mt(2).TTT, ...
                     mt(1).DD, mt(1).ZZ, mt(1).VVall, A0, P0);
else

  mt(1).zend = A0;
  mt(1).Pend = P0;

end


%% MAIN SAMPLE

% Loop over normal period and extra periods (if any)
for p = 2:pd
  if p > 2
    [zprev,Pprev] = augmentStates(mspec, mt(p).nant{:}, size(mt(p).TTT,2), mt(p-1).zend, mt(p-1).Pend);
  else
    zprev = mt(p-1).zend;
    Pprev = mt(p-1).Pend;
  end
  [mt(p).pyt, mt(p).zend, mt(p).Pend] = ...
    kalcvf2NaN( (mt(p).YY)', 1, zeros(size(mt(p).TTT,2),1), mt(p).TTT, ...
      mt(p).DD, mt(p).ZZ, mt(p).VVall, zprev, Pprev);
end

pyt = sum([mt(2:end).pyt]);



