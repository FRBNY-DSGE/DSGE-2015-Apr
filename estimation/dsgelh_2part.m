function pyt = dsgelh_2part(para,YY,nobs,nlags,nvar,mspec,npara,coint,cointadd,nant, antlags)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
%% transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[TTT,RRR,CCC,valid] = dsgesolv(mspec,para, nant);

iSplit = getState(mspec,0,'n_end+n_exo');
nstate = size(TTT,1);
nshocks = size(RRR,2)-nant;

% For 510,555, TTT/RRR/CCC for 510 are contained within those for 555
if any(mspec == [555 556 557 5571 558])
    TTT1 = zeros(nstate-(nant+1),nstate-(nant+1));  
    TTT1 = [TTT(1:iSplit,1:iSplit), ...
             TTT(1:iSplit,iSplit+(nant+1)+1:end);...
             TTT(iSplit+(nant+1)+1:end,1:iSplit),...
             TTT(iSplit+(nant+1)+1:end,iSplit+(nant+1)+1:end)];
    
    RRR1 = zeros(nstate-(nant+1),nshocks);
    RRR1 = [RRR(1:iSplit,1:nshocks); RRR(iSplit+(nant+1)+1:end,1:nshocks)];
    
    CCC1 = zeros(nstate-(nant+1),1);
    CCC1 = [CCC(1:iSplit,1); CCC(iSplit+(nant+1)+1:end,1)];
else 
    TTT1 = TTT;
    RRR1 = RRR;
    CCC1 = CCC;
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
%keyboard;
eval(['[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT1,RRR1,valid,para,nvar-nant,nlags,mspec,npara,coint,cointadd);']);

if retcode == 0
    % invalid parameterization
    pyt = -1E10;
    return;
end;

nstate  = size(TTT1,1);
nshocks = size(RRR1,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: compute log-likelihood using Kalman filter - written by Iskander
%%         note that Iskander's program assumes a transition equation written as:
%%         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t 
%%         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
%%         and  VV2 = cov(eps2_t,u_u) = RRR*VV
%%         define VVall as the joint variance of the two shocks VVall = var([eps2_t;u_t])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% truncate system matrices - these should NOT include the expectations

%taking the expectations out of the data
%nvar = nvar-nant;
YY1= YY(:, 1:nvar-nant);

HH = EE+MM*QQ*MM';
VV = QQ*MM';
VVall = [[RRR1*QQ*RRR1',RRR1*VV];[VV'*RRR1',HH]];

if any(CCC1 ~= 0)
  DDadd = ZZ*((eye(size(TTT1))-TTT1)\CCC1);
else
  DDadd = 0;
end

DD= DD+DDadd;


%% Define the initial mean and variance for the state vector
% (deals with nonstationarity on model by model basis)
[A0,P0] = lyap_nonstationary(mspec,para,TTT1,RRR1,QQ);

  zend = A0;
  Pend = P0;


%% Running Kalman Filter normal model (main sample)

% First, filter using normal model up to period T-antlags-1.
[pyt1,zend,Pend] = kalcvf2NaN(YY1(1:end-antlags-1,:)',1,zeros(nstate,1),TTT1,DD,ZZ,VVall,zend,Pend);

% Now change models to incorporate anticipated shocks.
eval(['[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT,RRR,valid,para,nvar,nlags,mspec,npara,coint,cointadd, nant);']);

HH = EE+MM*QQ*MM';
VV = QQ*MM';
VVall = [[RRR*QQ*RRR',RRR*VV];[VV'*RRR',HH]];

nstate = size(TTT,1);
nshocks = size(RRR,2);

[zprev,pprev] = augmentStates(mspec,nant,nstate,zend,Pend);
% Next, filter using ZB model for periods T-antlags:T
% These new filtered states include anticipated shocks and
% are consistent with the zero bound.
[pyt2,zend,pend] = kalcvf2NaN(YY(end-antlags:end,:)',1,zeros(nstate,1),TTT,DD,ZZ,VVall,zprev,pprev);

pyt = pyt1+pyt2;


