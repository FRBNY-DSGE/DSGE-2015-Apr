function [alpha_til,eta_til] = drawstates_dk02(YY,YY0,params,mspec,A0,P0,nant,antlags,nConditionalQuarters,zerobound,P0_ZB,ind_r, r_tl1, r_tlx)


% DRAWSTATES_DK02.M
% Originally written by DG.
% Last updated by SS 05/27/2016. This program was substantially modified to
% properly handle the zerobound period. Prior to this fix, the Kalman
% filtering that takes place below was not split into two regimes (pre-ZLB
% and ZLB) like we normally do.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is a simulation smoother based on Durbin and Koopman's
% "A Simple and Efficient Simulation Smoother for State Space Time Series
% Analysis" (Biometrika, 2002). The algorithm has been simplified for the
% case in which there is no measurement error, and the model matrices do
% not vary with time.Note that the state outputted by this program should
% equal the Kalman filtered state in expectation.

% Unlike other simulation smoothers (for example, that of Carter and Kohn,
% 1994), this method does not require separate draws for each period, draws
% of the state vectors, or even draws from a conditional distribution.
% Instead, vectors of shocks are drawn from the unconditional distribution
% of shocks, which is then corrected (via a Kalman Smoothing step), to
% yield a draw of shocks conditional on the data. This is then used to
% generate a draw of states conditional on the data. Drawing the states in
% this way is much more efficient than other methods, as it avoids the need
% for multiple draws of state vectors (requiring singular value
% decompositions), as well as inverting state covariance matrices
% (requiring the use of the computationally intensive and relatively
% erratic Moore-Penrose pseudoinverse).

% The state space is assumed to take the form:
% y(t) = ZZ*alpha(t) + DD
% alpha(t+1) = TTT*alpha(t) + RRR*eta(t+1)

% INPUTS:
% YY, the (Ny x Nt) matrix of main-sample observable data.
% YY0, the (Ny x Nt) matrix of pre-sample observable data.
% params, the vector of parameters
% mspec, the mspec number for the model of interest, necessary to solve
%       model
% A0, the (Nz x 1) time-invariant state vector.
% P0, the (Nz x Nz) time-invariant state covariance matrix.
% nant, a scalar for the zero bound specification indicating the
%       number of periods ahead the interest rate is fixed. When
%       zerobound=0, nant=0.
% antlags, a scalar for the zero bound specification indicating
%       the number of periods for which interest rate expectations have
%       been fixed
% nConditionalQuarters, a scalar for the number of periods of conditional
%       data included in YY matrix.
% P0_ZB is (Nz+nant x Nz+nant) time-invariant state vector for augmented
%       states (i.e., including the anticipated shock states). This is an
%       optional argument and is only necessary if you run with
%       zerobound=1.
% ind_r, r_tl1, r_tlx are state indices. These are optional arguments and
%       only necessary if you run with zerobound=1.

% OUTPUTS:
% alpha_til, the vector of smoothed states
% eta_til, the vector of smoothed shocks


%% Get system matrices
% Note: If zerobound==1, here we get the augmented system matrices.
[TTT,RRR,CCC,valid_a] = dsgesolv(params,mspec,nant);
[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = getmeasur(mspec,TTT,RRR,valid_a,params,...
    size(YY,1),0,size(A0,1),0,0,nant);

%% Get matrix dimensions
Ny  = size(YY,1);  % number of observables
Nt0 = size(YY0,2); % number of periods in the pre-sample
Nt  = size(YY,2);  % number of periods in the main sample
Nz  = size(TTT,1); % number of states
Ne  = size(RRR,2); % number of shocks

%% Produce "fake" states and observables (a+ and y+)

% initialize matrices
alpha_all_plus = nan(Nz,Nt0+Nt);
YY_all_plus = nan(Ny,Nt0+Nt);

% Draw initial state a_0+
if zerobound
    [U,D,V] = svd(P0_ZB);
else
    [U,D,V] = svd(P0);
end
ap_t = U*sqrt(D)*randn(Nz,1);

% Draw a sequence of shocks eta+
eta_all_plus = sqrt(QQ)*randn(Ne,Nt0+Nt);

% Set nant shocks to 0 for non-ZLB time periods
if nant > 0
    eta_all_plus(end-nant+1:end,1:end-antlags-(nConditionalQuarters+1)) = 0;
end

% Iterate forward state and observation equations
for t = 1:Nt0+Nt
    ap_t = TTT*ap_t + RRR*eta_all_plus(:,t);
    alpha_all_plus(:,t) = ap_t;
    YY_all_plus(:,t) = ZZ*ap_t+DD;
end

% replace fake data with NaNs wherever actual data has NaNs
YY_all= [YY0 YY];
YY_all_plus(isnan(YY_all))=nan;

%% Compute y* = y - y+ - D
YY_star = YY_all - YY_all_plus;

%% Filter and smooth

if zerobound
    %%%%% Step 1: Kalman filter under pre-ZLB regime  %%%%%
    
    % get pre-ZLB system matrices
    [TTT_preZLB,RRR_preZLB,CCC_preZLB,valid_a] = dsgesolv(params,mspec);
    Nz_preZLB  = size(TTT_preZLB,1);
    [ZZ_preZLB,DD_preZLB,DDcointadd,QQ_preZLB,EE_preZLB,MM_preZLB,retcode] = ...
        getmeasur(mspec,TTT_preZLB,RRR_preZLB,valid_a,params,Ny-nant,0,Nz_preZLB,0,0);
    
    % get time-invariant variance matrix of shocks (NB: we assume no
    % measurement error here)
    var_preZLB = zeros(Ny-nant+Nz_preZLB,Ny-nant+Nz_preZLB);
    var_preZLB(1:Nz_preZLB,1:Nz_preZLB) = RRR_preZLB*QQ_preZLB*RRR_preZLB';
    
    % Kalman filter over pre-ZLB era
    [L,zend,pend,pred,vpred] = kalcvf2NaN(YY_star(1:Ny-nant,1:end-antlags-1-nConditionalQuarters),0,...
        zeros(size(CCC_preZLB)),TTT_preZLB,zeros(size(DD_preZLB)),ZZ_preZLB,...
        var_preZLB,zeros(size(A0)),P0);
    
    %%%%% Step 2: Kalman filter under ZLB regime  %%%%%%
    
    % augment states to include anticipated policy shocks
    ExpFFR=YY_star(end-nant+1:end,end-antlags-nConditionalQuarters:end)';
    if is2part(mspec)
        [zprev,pprev,MM_ant,EE_ant,ZZ_e,DD_e,YY_ant] = augmentFilter(r_tl1,r_tlx,nant,...
            antlags+nConditionalQuarters,ind_r,zend,pend,MM(1:end-nant,:),EE(1:end-nant,1:end-nant),...
            ZZ(1:end-nant,:),zeros(size(DD(1:end-nant))),YY_star(1:end-nant,:)',...
            TTT,ExpFFR);
    else
        [zprev,pprev,MM_ant,EE_ant,ZZ_e,DD_e,YY_ant] = augmentFilter(r_tl1,r_tlx,nant,...
            antlags,ind_r,zend,pend,MM,EE,ZZ,zeros(size(DD)),YY_star,TTT,ExpFFR);
    end
    var = zeros(Ny+Nz,Ny+Nz);
    var(1:Nz,1:Nz) = RRR*QQ*RRR';
    
    % Kalman filter over ZLB era
    [L_ant,zend,pend,pred_ant,vpred_ant] = kalcvf2NaN(YY_ant',0,zeros(size(CCC)),TTT,...
        zeros(size(DD_e)),ZZ_e,var,zprev,pprev);
    
    %%%%% Step 3: Kalman smooth over everything %%%%%
    [pred_all,vpred_all,YY_all,A0_,P0_] = augmentSmoother(Nz,r_tlx,r_tl1,nant,antlags+nConditionalQuarters,...
        [],pred,pred_ant,[],vpred,vpred_ant,[],YY_star',YY_ant,A0,P0);
    
    [alpha_hat_star,eta_hat_star] = kalsmth_k93(A0_,P0_,YY_all',pred_all,vpred_all,...
        TTT,RRR,QQ,ZZ_e,zeros(size(DD_e)),nant,antlags,nConditionalQuarters>0,nConditionalQuarters);
    
else
    % get time-invariant variance matrix of shocks (NB: we assume no
    % measurement error here)
    var = zeros(Ny+Nz,Ny+Nz);
    var(1:Nz,1:Nz) = RRR*QQ*RRR';
    
    % Kalman filter
    [L,zend,pend,pred,vpred] = kalcvf2NaN(YY_star,0,zeros(size(CCC)),TTT,zeros(size(DD)),...
        ZZ,var,A0,P0);
    
    % Kalman smooth
    [alpha_hat_star,eta_hat_star] = kalsmth_k93(A0,P0,YY_star,pred,vpred,TTT,RRR,QQ,...
        ZZ,zeros(size(DD)),nant,antlags,nConditionalQuarters>0,nConditionalQuarters);
end

%% Compute draw (states and shocks)
alpha_til = alpha_all_plus + alpha_hat_star;
eta_til = eta_all_plus + eta_hat_star;



