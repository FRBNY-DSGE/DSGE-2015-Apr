function [para,para_names,para_mask,para_fix,npara,polipar,polivalue,bounds] = mspec_parameters_990(subspec,dataset)

%% names of parameters
% the parameters in para_names must be lsited in the order they appear
% below and in the model's priorfile

nantpad = 20;

para_names = [...
  { '\alpha';'\zeta_p';'\iota_p';'\Upsilon';'\Phi';'S''''';'h';'psi';'\nu_l';'\zeta_w';...
    '\iota_w';'\beta';...
    '\psi_1';'\psi_2';'\psi_3';'\pi^*';'\sigma_{c}';'\rho';'F(\omega)';'spr_*';'\zeta_{sp}';'\gamma_*';...
    '\gamma';'Lmean';...
    ...
    '\rho_{g}';'\rho_{b}';'\rho_{\mu}';'\rho_{z}';'\rho_{\lambda_f}';'\rho_{\lambda_w}';...
    '\rho_{rm}';'\rho_{sigw}';'\rho_{mue}';'\rho_{gamm}';'\rho_{pi}^*';'\rho_{lr}';'\rho_{z^p}';'\rho_{tfp}';'\rho_{gdpdef}';'\rho_{pce}';...
    ...
    '\sigma_{g}';'\sigma_{b}';'\sigma_{\mu}';'\sigma_{z}';'\sigma_{\lambda_f}';'\sigma_{\lambda_w}';...
    '\sigma_{rm}';'\sigma_{sigw}';'\sigma_{mue}';'\sigma_{gamm}';'\sigma_{pi}^*';'\sigma_{lr}';'\sigma_{z^p}';'\sigma_{tfp}';'\sigma_{gdpdef}';'\sigma_{pce}';...
  }; ...
    ...
    arrayfun(@(i_) sprintf('\\sigma_{ant%d}', i_), [1:nantpad]', 'un', 0);
    ...
    {'\eta_{gz}';'\eta_{\lambda_f}';'\eta_{\lambda_w}'; 'i_{\alpha}^{model}'; '\Gamma_{gdpdef}'; '\delta_{gdpdef}'} ...
  ];

%% values for parameters

% fixing the parameters

para        = zeros(100,1);
para_mask   = zeros(100,1);   %npara x 1, set corresponding element to 1 to fix the parameters
para_fix    = zeros(100,1);   %npara x 1, contains values for fixed parameters

if subspec == 93
  para(1)     = 0.3000;    %para_mask(1) = 1;   para_fix(1) = para(1);              %% alp;         1
else
  para(1)     = 0.1596;    %para_mask(1) = 1;   para_fix(1) = para(1);              %% alp;         1
end
para(2)     = 0.8940;  %para_mask(2) = 1;   para_fix(2) = para(2);              %% zeta_p;      2
para(3)     = 0.1865;  %para_mask(3) = 1;   para_fix(3) = para(3);              %% iota_p;      3  
para(4)     = 1.0000;  para_mask(4) = 1;   para_fix(4) = para(4);              %% ups; .1      5
para(5)     = 1.1066; %para_mask(5) = 1;   para_fix(5) = para(5);              %% Bigphi;      6
para(6)     = 2.7314; %para_mask(6) = 1;   para_fix(6) = para(6);              %% s2;          7
para(7)     = 0.5347; %para_mask(7) = 1;   para_fix(7) = para(7);              %% h;           8
para(8)     = 0.6862; %para_mask(8) = 1;   para_fix(8) = para(8);              %% ppsi;        9
para(9)     = 2.5975; %para_mask(9) = 1;   para_fix(9) = para(9);              %% nu_l         10
para(10)    = 0.9291; %para_mask(10) = 1;  para_fix(10) = para(10);            %% zeta_w;      11
para(11)    = 0.2992; %para_mask(11) = 1;  para_fix(11) = para(11);            %% iota_w;      12
para(12)    = 0.1402; %para_mask(12) = 1;  para_fix(12) = para(12);            %% bet;         14
para(13)    = 1.3679; %para_mask(13) = 1;  para_fix(13) = para(13);            %% psi1;        15
para(14)    = 0.0388; %para_mask(14) = 1;  para_fix(14) = para(14); %% psi2;        16
para(15)    = 0.2464; %para_mask(15) = 1;  para_fix(15) = para(15); %% psi3;        17
if subspec == 2
  para(16)    = 0.5000; para_mask(16) = 1;  para_fix(16) = para(16); %% pistar;      18
else
  para(16)    = 1.1121; %para_mask(16) = 1;  para_fix(16) = para(16); %% pistar;      18
end
para(17)    = 0.8719; %para_mask(17) = 1;  para_fix(17) = para(17);     %% sigmac;      19
para(18)    = 0.7126; %para_mask(18) = 1;  para_fix(18) = para(18);     %% rho;      19

%Financial Frictions Parameters
para(19)    = 0.0300;  para_mask(19) = 1;  para_fix(19) = para(19);              %% Fom
para(20)    = 1.7444;  %para_mask(20) = 1;  para_fix(20) = para(20);              %% st st spread
para(21)    = 0.0559;  %para_mask(21) = 1;  para_fix(21) = para(21);              %% zeta_spb0.050;
para(22)    = 0.9900;  para_mask(22) = 1;  para_fix(22) = para(22);              %% gammstar

npara       = 22;           %number of parameters

%% exogenous processes - level

para(npara+1) = 0.3673; %para_mask(npara+1) = 1;    para_fix(npara+1) = para(npara+1);      %% gam;     
para(npara+2) = -45.9364;    %para_mask(npara+2) = 1;    para_fix(npara+2) = para(npara+2);      %% Lmean;

npara = npara+2;
    
% exogenous processes - autocorrelation
        
para(npara+1)   = 0.9863; %para_mask(npara+1) = 1;       para_fix(npara+1) = para(npara+1);      %% rho_g
para(npara+2)   = 0.9410; %para_mask(npara+2) = 1;       para_fix(npara+2) = para(npara+2);      %% rho_b;
para(npara+3)   = 0.8735; %para_mask(npara+3) = 1;       para_fix(npara+3) = para(npara+3);      %% rho_mu
para(npara+4)   = 0.9446; %para_mask(npara+4) = 1;       para_fix(npara+4) = para(npara+4);      %% rho_z;
para(npara+5)   = 0.8827; %para_mask(npara+5) = 1;       para_fix(npara+5) = para(npara+5);      %% rho_laf;
para(npara+6)   = 0.3884; %para_mask(npara+6) = 1;       para_fix(npara+6) = para(npara+6);      %% rho_law;
para(npara+7)   = 0.2135; %para_mask(npara+7) = 1;       para_fix(npara+7) = para(npara+7);      %% rho_rm;
para(npara+8)   = 0.9898; %para_mask(npara+8) = 1;       para_fix(npara+8) = para(npara+8);      %% rho_sigw
para(npara+9)   = 0.7500;  para_mask(npara+9) = 1;       para_fix(npara+9) = para(npara+9);      %% rho_mue
para(npara+10)  = 0.7500;  para_mask(npara+10) = 1;      para_fix(npara+10) = para(npara+10);    %% rho_gamm
para(npara+11)  = 0.9900;  para_mask(npara+11) = 1;       para_fix(npara+11) = para(npara+11);      %% rho_pist;
para(npara+12)  = 0.6936; %para_mask(npara+12) = 1;       para_fix(npara+12) = para(npara+12);      %% rho_lr;
para(npara+13)  = 0.8910; %para_mask(npara+13) = 1;       para_fix(npara+13) = para(npara+13);      %% rho_zp;
para(npara+14)  = 0.1953; %para_mask(npara+14) = 1;       para_fix(npara+14) = para(npara+14);      %% rho_tfp;
if subspec == 94  
  para(npara+15)  = 0.0000; para_mask(npara+15) = 1;       para_fix(npara+15) = para(npara+15);      %% rho_gdpdef;
  para(npara+16)  = 0.0000; para_mask(npara+16) = 1;       para_fix(npara+16) = para(npara+16);      %% rho_pce;
else
  para(npara+15)  = 0.5379; %para_mask(npara+15) = 1;       para_fix(npara+15) = para(npara+15);      %% rho_gdpdef;
  para(npara+16)  = 0.2320; %para_mask(npara+16) = 1;       para_fix(npara+16) = para(npara+16);      %% rho_pce;
end

npara = npara+16;
               
%% exogenous processes - standard deviation    
para(npara+1)  = 2.5230; %para_mask(npara+1) = 1;       para_fix(npara+1) = para(npara+1);      %% sig_g;
para(npara+2)  = 0.0292; %para_mask(npara+2) = 1;       para_fix(npara+2) = para(npara+2);      %% sig_b;
para(npara+3)  = 0.4559;  %para_mask(npara+3) = 1;       para_fix(npara+3) = para(npara+3);      %% sig_mu;
para(npara+4)  = 0.6742;  %para_mask(npara+4) = 1;       para_fix(npara+4) = para(npara+4);      %% sig_z;
para(npara+5)  = 0.1314;  %para_mask(npara+5) = 1;       para_fix(npara+5) = para(npara+5);      %% sig_laf;
para(npara+6)  = 0.3864;  %para_mask(npara+6) = 1;       para_fix(npara+6) = para(npara+6);      %% sig_law;
para(npara+7)  = 0.2380;  %para_mask(npara+7) = 1;       para_fix(npara+7) = para(npara+7);      %% sig_rm;
para(npara+8)  = 0.0428;    %para_mask(npara+8) = 1;      para_fix(npara+8) = para(npara+8);     %% sig_sigw
para(npara+9)  = 0     ;   para_mask(npara+9) = 1;      para_fix(npara+9) = para(npara+9);      %% sig_mue
para(npara+10) = 0     ;   para_mask(npara+10) = 1;     para_fix(npara+10) = para(npara+10);      %% sig_gamm
para(npara+11) = 0.0269;  %para_mask(npara+11) = 1;       para_fix(npara+11) = para(npara+11);      %% sig_pist;
para(npara+12) = 0.1766;  %para_mask(npara+12) = 1;       para_fix(npara+12) = para(npara+12);      %% sig_lr;
para(npara+13) = 0.1662;  %para_mask(npara+13) = 1;       para_fix(npara+13) = para(npara+13);      %% sig_zp;
if subspec == 9 
  para(npara+14)  = 0.8000;  para_mask(npara+14) = 1;       para_fix(npara+14) = para(npara+14);      %% sig_tfp;
elseif any(subspec == [2 90 93 94])
  para(npara+14)  = 0.9391;  %para_mask(npara+14) = 1;       para_fix(npara+14) = para(npara+14);      %% sig_tfp;
else
  error('Not implemented.');
end
para(npara+15)  = 0.1575;  %para_mask(npara+15) = 1;       para_fix(npara+15) = para(npara+15);      %% sig_gdpdef;
para(npara+16)  = 0.0999;  %para_mask(npara+16) = 1;       para_fix(npara+16) = para(npara+16);      %% sig_pce;

npara = npara+16;


% Standard Deviations of the anticipated policy shocks
for i = 1:nantpad
  eval(strcat('para(npara +',num2str(i),') = 0.20;'));
  if i >=13
    eval(strcat('para_mask(npara +',num2str(i),') = 1;'));
    eval(strcat('para_fix(npara +',num2str(i),') = 0;'));
  end
end
npara = npara+nantpad;

para(npara+1) = 0.8400;  %para_mask(npara+12) = 1;      para_fix(npara+12) = para(npara+12);      %% eta_gz;
para(npara+2) = 0.7892;  %para_mask(npara+13) = 1;      para_fix(npara+13) = para(npara+13);      %% eta_laf;
para(npara+3) = 0.4226;  %para_mask(npara+14) = 1;      para_fix(npara+14) = para(npara+14);    %% eta_law;
npara = npara+3;

if subspec == 9 
  para(npara+1) = 1.0000;  para_mask(npara+1) = 1;      para_fix(npara+1) = para(npara+1);    %% modelalp_ind, datasets 712;
elseif any(subspec == [2 90 93 94])
  para(npara+1) = 0.0000;  para_mask(npara+1) = 1;      para_fix(npara+1) = para(npara+1);    %% modelalp_ind datasets 714;
else
  error('Not implemented.');
end
para(npara+2) = 1.0354;       %para_mask(npara+3) = 3;      para_fix(npara+3) = para(npara+3);    %% gamm_gdpdef
para(npara+3) = 0.0181;       %para_mask(npara+4) = 4;      para_fix(npara+4) = para(npara+4);    %% del_gdpdef
npara = npara+3;

para = para(1:npara);
para_mask = para_mask(1:npara);
para_fix = para_fix(1:npara);

%% initialize parameters at fixed value or 0   

para = para.*(1-para_mask)+para_fix.*para_mask;

%% identify policy and non policy parameters - I don't think we use these

polipar = 16;
polivalue = 2;

%% bounds for MH

bounds = zeros(100,2);

bounds(1,:) = [1E-5 .999];      %% alp;         1
bounds(2,:) = [1E-5 .999];      %% zeta_p;      2
bounds(3,:) = [1E-5 .999];      %% iota_p;      3
bounds(4,:) = [0    10];        %% ups;         5
bounds(5,:) = [1    10];        %% Bigphi;      6
bounds(6,:) = [-15   15];       %% s2;          7
bounds(7,:) = [1E-5 .999];      %% h;           8
bounds(8,:) = [1E-5 .999];      %% ppsi;        9
bounds(9,:) = [1E-5 10];        %% nu_l         10
bounds(10,:) = [1E-5 .999];     %% zeta_w;      11
bounds(11,:) = [1E-5 .999];     %% iota_w;      12
bounds(12,:) = [1E-5    10];    %% bet;         14
bounds(13,:) = [1E-5    10];    %% psi1;        15
bounds(14,:) = [-.5     .5];    %% psi2;        16
bounds(15,:) = [-.5     .5];    %% psi3;        17
bounds(16,:) = [1E-5  10];      %% pistar;      18
bounds(17,:) = [1E-5    10];    %% sigmac;      19
bounds(18,:) = [1E-5    .999];  %% rho;      19

bounds(19,:) = [1E-5   .99999]; %% Fom
bounds(20,:) = [0 100]; %% st st spread
bounds(21,:) = [1E-5   .99999]; %% zeta_sp
bounds(22,:) = [1E-5   .99999]; %% gammstar

npara = 22;

%% exogenous processes - level

bounds(npara+1,:) = [-5   5];        %% gam;     22
bounds(npara+2,:) = [-1000   1000];  %% Lstar;   23

npara = npara+2;    

%% exogenous processes - autocorrelation

bounds(npara+1,:) = [1E-5   .999];    %% rho_g;
bounds(npara+2,:) = [1E-5   .999];    %% rho_b;
bounds(npara+3,:) = [1E-5   .999];    %% rho_mu;
bounds(npara+4,:) = [1E-5   .999];    %% rho_z;
bounds(npara+5,:) = [1E-5   .999];    %% rho_laf;
bounds(npara+6,:) = [1E-5   .999];    %% rho_law;
bounds(npara+7,:) = [1E-5   .999];    %% rho_rm;

bounds(npara+8,:) = [1E-5   .99999]; %% rho_sigw
bounds(npara+9,:) = [1E-5   .99999]; %% rho_mue
bounds(npara+10,:)= [1E-5   .99999]; %% rho_gamm

bounds(npara+11,:) = [1E-5   .999];    %% rho_pist;
bounds(npara+12,:) = [1E-5   .999];    %% rho_lr;

bounds(npara+13,:) = [1E-5   .999];    %% rho_zp;
bounds(npara+14,:) = [1E-5   .999];    %% rho_tfp;
bounds(npara+15,:) = [1E-5   .999];    %% rho_gdpdef;
bounds(npara+16,:) = [1E-5   .999];    %% rho_pce;
npara = npara+16;

%% exogenous processes - standard deviation

bounds(npara+1,:) = [1E-8     5];       %% sig_g (sig_zP if mspec == 8);
bounds(npara+2,:) = [1E-8     5];       %% sig_b;
bounds(npara+3,:) = [1E-8     5];       %% sig_mu;
bounds(npara+4,:) = [1E-8     5];       %% sig_z;
bounds(npara+5,:) = [1E-8     5];       %% sig_laf;
bounds(npara+6,:) = [1E-8     5];       %% sig_law;
bounds(npara+7,:) = [1E-8     5];       %% sig_rm;

bounds(npara+8,:) = [1E-7   100];       %% sig_sigw
bounds(npara+9,:) = [1E-7   100];       %% sig_mue
bounds(npara+10,:) = [1E-7   100];      %% sig_gamm

bounds(npara+11,:) = [1E-8     5];       %% sig_pist;
bounds(npara+12,:) = [1E-8     10];       %% sig_lr;

bounds(npara+13,:) = [1E-8     5];       %% sig_z;
bounds(npara+14,:) = [1E-8     5];       %% sig_tfp;
bounds(npara+15,:) = [1E-8     5];       %% sig_gdpdef;
bounds(npara+16,:) = [1E-8     5];       %% sig_pce;
npara = npara+16;


%% Standard Deviations of the Anticipated Shocks
for i = 1:nantpad
    eval(strcat('bounds(npara +',num2str(i),',:) = [1E-7   100];'));
end
npara = npara+nantpad;

bounds(npara+1,:) = [1E-5  .999];      %% eta_gz;
bounds(npara+2,:) = [1E-5  .999];      %% eta_laf;
bounds(npara+3,:) = [1E-5  .999];      %% eta_law;
npara = npara+3;

bounds(npara+1,:) = [0.000 1.00];      %% modelalp_ind;
bounds(npara+2,:) = [-10 10];          %% gamm_gdpdef
bounds(npara+3,:) = [-9.1 9.1];        %% del_gdpdef

npara = npara+3;


bounds = bounds(1:npara,:);














