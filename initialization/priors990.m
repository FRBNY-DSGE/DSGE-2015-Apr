% Define Prior parameters
% pshape is 1: BETA(mean,stdd)
%           2: GAMMA(mean,stdd)
%           3: NORMAL(mean,stdd)
%           4: INVGAMMA(s^2,nu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loose lambda_f prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prior = pri990

prior = zeros(100,3);
nantpad = 20;

prior(1,:) = [.30   .05  3]; %% alp - beta -changed
prior(2,:) = [.5    .1  1];  %% zeta_p - beta -changed
prior(3,:) = [.5    .15 1];  %% iota_p - betaa
prior(4,:) = [1     .5 2];   %% ups - gamma
prior(5,:) = [1.25  .12 3];  %% Bigphi - gamma
prior(6,:) = [4    1.5  3];  %% s2 - gamma
prior(7,:) = [.7    .1  1];  %% h - beta -changed
prior(8,:) = [.5    .15 1];  %% ppsi - gamma
prior(9,:) = [2    .75 3];  %% nu_l - gamma - new mean
prior(10,:) = [.5   .1  1];  %% zeta_w - beta -changed
prior(11,:) = [.5   .15  1]; %% iota_w - beta
prior(12,:) = [.25   .1  2]; %% bet - beta
prior(13,:) = [1.5  .25  3];  %% psi1 - gamma -changed
prior(14,:) = [.12   .05  3];  %% psi2 - gamma
prior(15,:) = [.12   .05  3];  %% psi3 - gamma
prior(16,:) = [.75  .4  2];  %% pistar - gamma
prior(17,:) = [1.5  0.37  3];  %% sigmac - normal
prior(18,:) = [.75  0.10  1];  %% rho

prior(19,:) = [.03 .01  1];  %% Fom - beta
prior(20,:) = [2   .1   2];  %% st st spread - gamma
prior(21,:) = [.05 .005  1];  %% zeta_sp - beta
prior(22,:) = [.99 .002 1];  %% gammstar - beta


npara = 22; 

%% exogenous processes - level
prior(npara+1,:) = [.4 .1 3];      %% gam - normal
prior(npara+2,:) = [-45  5  3];      %% Lmean - normal
%prior(npara+3,:) = [0.71 1  3];      %% zconst - normal

npara = npara+2;    

%% exogenous processes - autocorrelation -CHANGE TO STANDARD!
prior(npara+1,:) = [.5   .2  1];    %% rho_g - beta
prior(npara+2,:) = [.5   .2  1];    %% rho_b - beta
prior(npara+3,:) = [.5   .2  1];    %% rho_mu - beta
prior(npara+4,:) = [.5   .2  1];    %% rho_z - beta
prior(npara+5,:) = [.5   .2  1];    %% rho_laf - beta
prior(npara+6,:) = [.5   .2  1];    %% rho_law - beta
prior(npara+7,:) = [.5   .2  1];    %% rho_rm - beta

prior(npara+8,:) = [.75  .15  1];   %% rho_sigw - beta
prior(npara+9,:) = [.75  .15  1];   %% rho_mue - beta
prior(npara+10,:) = [.75  .15  1];   %% rho_gamm - beta
prior(npara+11,:) = [.5   .2  1];    %% rho_pist - beta
prior(npara+12,:) = [.5   .2  1];    %% rho_lr - beta
prior(npara+13,:) = [.5   .2  1];    %% rho_zp - beta
prior(npara+14,:) = [.5   .2  1];    %% rho_tfp - beta
prior(npara+15,:) = [.5   .2  1];    %% rho_gdpdef - beta
prior(npara+16,:) = [.5   .2  1];    %% rho_pce - beta

npara = npara+16;

%% exogenous processes - standard deviation
prior(npara+1,:) = [0.10, 2.00, 4]; %% sig_g;
prior(npara+2,:) = [0.10, 2.00, 4]; %% sig_b;
prior(npara+3,:) = [0.10, 2.00, 4]; %% sig_mu;
prior(npara+4,:) = [0.10, 2.00, 4]; %% sig_z;
prior(npara+5,:) = [0.10, 2.00, 4]; %% sig_laf;
prior(npara+6,:) = [0.10, 2.00, 4]; %% sig_law;
prior(npara+7,:) = [0.10, 2.00, 4]; %% sig_rm;

prior(npara+8,:) = [.20/4  4.00    4]; %% sig_sigw;
prior(npara+9,:) = [.20/4  4.00    4]; %% sig_mue;
prior(npara+10,:) = [.01  4.00    4]; %% sig_gamm;

prior(npara+11,:) = [0.03, 6, 4]; %% sig_pist;
prior(npara+12,:) = [0.75, 2, 4]; %% sig_lr;

prior(npara+13,:) = [0.10, 2.00, 4]; %% sig_zp;
prior(npara+14,:) = [0.10, 2.00, 4]; %% sig_tfp;
prior(npara+15,:) = [0.10, 2.00, 4]; %% sig_gdpdef;
prior(npara+16,:) = [0.10, 2.00, 4]; %% sig_pce;

npara = npara+16;

for i = 1:nantpad
    eval(strcat('prior(npara +',num2str(i),',:) = [.2  4.00    4];'));
end
npara = npara+nantpad;

prior(npara+1,:) = [0.50, 0.20, 1]; %% eta_gz;
prior(npara+2,:) = [0.50, 0.20, 1]; %% eta_laf;
prior(npara+3,:) = [0.50, 0.20, 1]; %% eta_law;

npara = npara+3;

prior(npara+1,:) = [0.50, 0.20, 1]; %% modelalp_ind;
prior(npara+2,:) = [1.00, 2,    3]; %% gamm_gdpdef;
prior(npara+3,:) = [0.00, 2,    3]; %% del_gdpdef;
npara = npara+3;

prior = prior(1:npara,:);
