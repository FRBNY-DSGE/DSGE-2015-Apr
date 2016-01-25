function [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur990(TTT,RRR,valid,para,nvar,nlags,mspec,npara,coint,cointadd,nant,varargin);
%% description:
%% solution to DSGE model - delivers transition equation for the state variables  S_t
%% transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
%% define the measurement equation: X_t = ZZ S_t +D+u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

if length(varargin) > 0
    subspec = varargin{1};
end

retcode = 1;

if valid < 1;
    retcode = 0;

    ZZ = [];
    DD = [];
    QQ = [];
    EE = [];
    MM = [];    
    DDcointadd = [];
    
    return
end

nstate = size(TTT,1);
DDcointadd = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: assign names to the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
getPara_script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: assign names to the columns of GAM0, GAM1 -- state variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(strcat('states',num2str(mspec)));

%% additional states
y_t1 = n_end+n_exo+n_exp+1;
c_t1 = n_end+n_exo+n_exp+2;
i_t1 = n_end+n_exo+n_exp+3;
w_t1 = n_end+n_exo+n_exp+4;
pi_t1 = n_end+n_exo+n_exp+5;     % Add lagged mc_t state
L_t1  = n_end+n_exo+n_exp+6;     % Add lagged L_t state
Et_pi_t = n_end+n_exo+n_exp+7;   % Add forward looking expected infl.
lr_t = n_end+n_exo+n_exp+8;      % Add in measurement error for long run series
tfp_t = n_end+n_exo+n_exp+9;     % Add measurement error for Fernald TFP series
e_gdpdef = n_end+n_exo+n_exp+10; % Add measurement error for GDP Deflator
e_pce = n_end+n_exo+n_exp+11;    % Add measurement error for Core PCE
u_t1 = n_end+n_exo+n_exp+12;     % Add measurement error for Fernald TFP series

if nstate ~= (n_end+n_exo+n_exp+12)

    retcode = 0;

    yyyyd = zeros(nvar,nvar);
    xxyyd = zeros(1+nlags*nvar,nvar);
    xxxxd = zeros(1+nlags*nvar,1+nlags*nvar);

    disp('\n\n number of states does not match in vaprio\n');
    return
end

if ~exist('nant','var')
    nvar = 12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: assign measurement equation : X_t = ZZ*S_t + DD + u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create system matrices for state space model
ZZ = zeros(nvar+coint,nstate);
%% constant
DD = zeros(nvar+coint,1);
%% cov(eps_t,u_t) = VV
MM = zeros(nvar+coint,nex);   
%% var(eta_t) = EE
EE = zeros(nvar+coint,nvar+coint);
%% var(eps_t) = QQ
QQ =  zeros(nex,nex);



%% Output growth - Quarterly! 
ZZ(1,y_t) = 1;
ZZ(1,y_t1) = -1;
ZZ(1,z_t) = 1;
DD(1) = 100*(exp(zstar)-1);

%% Hoursg
ZZ(2,L_t) = 1;
DD(2) = Lmean;

%% Labor Share/real wage growth
%     ZZ(3,L_t) = 1;
%     ZZ(3,w_t) = 1;
%     ZZ(3,y_t) = -1;
%     DD(3) = 100*log((1-alp)/(1+laf));

ZZ(3,w_t) = 1;
ZZ(3,w_t1) = -1;
ZZ(3,z_t) = 1;

DD(3) = 100*(exp(zstar)-1);

%% Inflation (GDP Deflator)
ZZ(4,pi_t) = gamm_gdpdef;
ZZ(4,e_gdpdef) = 1;
DD(4) = 100*(pistar-1) + del_gdpdef;

%% Inflation (Core PCE)
ZZ(5,pi_t) = 1;
ZZ(5,e_pce) = 1;
DD(5) = 100*(pistar-1);

%% Nominal interest rate
ZZ(6,R_t) = 1;
DD(6) = Rstarn;

%% Consumption Growth
ZZ(7,c_t) = 1; 
ZZ(7,c_t1) = -1; 
ZZ(7,z_t) = 1; 
DD(7) = 100*(exp(zstar)-1);

%% Investment Growth
ZZ(8,i_t) = 1; 
ZZ(8,i_t1) = -1; 
ZZ(8,z_t) = 1; 
DD(8) = 100*(exp(zstar)-1);

%% Spreads
ZZ(9,E_Rktil) = 1;
ZZ(9,R_t) = -1;
DD(9) = 100*log(sprd);

%% 10 yrs infl exp

  % Create a TTT matrix that has the level rows/cols (i.e. the part with the
  % unit root) stripped out
  % strip_out = [zlev_t];
  % TTT_tmp = TTT; TTT_tmp(strip_out,:) = []; TTT_tmp(:,strip_out) = [];
  % replace_inds = reshape(1:nstate^2, nstate, nstate);
  % replace_inds(strip_out,:) = [];
  % replace_inds(:,strip_out) = [];
  % replace_inds = replace_inds(:);

  % Invert that matrix
  % TTT10_tmp = (1/40)*((eye(size(TTT_tmp,1)) - TTT_tmp)\(TTT_tmp - TTT_tmp^41));

  % Construct TTT10 matrix
  TTT10 = (1/40)*((eye(size(TTT,1)) - TTT)\(TTT - TTT^41));

ZZ(10,:) =  TTT10(pi_t,:);
DD(10) = 100*(pistar-1);

%% Long Rate
ZZ(11,:) = ZZ(6,:)*TTT10;
ZZ(11,lr_t) = 1;
DD(11) = Rstarn;

%% TFP
ZZ(12,z_t) = (1-alp)*modelalp_ind + 1*(1-modelalp_ind);
ZZ(12,tfp_t) = 1;
ZZ(12,u_t)  = alp/( (1-alp)*(1-modelalp_ind) + 1*modelalp_ind );
ZZ(12,u_t1) = -(alp/( (1-alp)*(1-modelalp_ind) + 1*modelalp_ind) );

QQ(g_sh,g_sh) = sig_g^2;
QQ(b_sh,b_sh) = sig_b^2;
QQ(mu_sh,mu_sh) = sig_mu^2;
QQ(z_sh,z_sh) = sig_z^2;
QQ(laf_sh,laf_sh) = sig_laf^2;
QQ(law_sh,law_sh) = sig_law^2;
QQ(rm_sh,rm_sh) = sig_rm^2;
QQ(sigw_sh,sigw_sh) = sig_sigw^2;  
QQ(mue_sh,mue_sh) = sig_mue^2;  
QQ(gamm_sh,gamm_sh) = sig_gamm^2;  
QQ(pist_sh,pist_sh) = sig_pist^2;
QQ(lr_sh,lr_sh) = sig_lr^2;
QQ(zp_sh,zp_sh) = sig_zp^2;
QQ(tfp_sh,tfp_sh)=sig_tfp^2;
QQ(gdpdef_sh,gdpdef_sh)=sig_gdpdef^2;
QQ(pce_sh,pce_sh)=sig_pce^2;

if exist('nant','var')
  if nant > 0
    % These lines set the standard deviations for the anticipated shocks. They
    % are here no longer calibrated to the std dev of contemporaneous shocks,
    % as we had in 904
    for i = 1:nant
        eval(strcat('QQ(rm_shl',num2str(i),',rm_shl',num2str(i),') = sig_rm_ant(', num2str(i), ')^2;'));
    end

    for i_ = 1:nant
      ZZ(12 + i_,:) = ZZ(6,:)*(TTT^i_);
      DD(12 + i_) = Rstarn;
    end
  end
end
%keyboard;
