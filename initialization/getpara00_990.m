function [alp,zeta_p,iota_p,del,ups,Bigphi,s2,h,ppsi,nu_l,zeta_w,iota_w,law,laf,bet,Rstarn,psi1,psi2,psi3,pistar,sigmac,rho,epsp,epsw...
    gam,Lmean,Lstar,gstar,rho_g,rho_b,rho_mu,rho_z,rho_laf,rho_law,rho_rm,rho_sigw,rho_mue,rho_gamm,rho_pist,rho_lr,rho_zp,rho_tfp,rho_gdpdef,rho_pce,...
    sig_g,sig_b,sig_mu,sig_z,sig_laf,sig_law,sig_rm,sig_sigw,sig_mue,sig_gamm,sig_pist,sig_rm_ant,sig_lr,sig_zp,sig_tfp,sig_gdpdef,sig_pce,...
    eta_gz,eta_laf,eta_law,modelalp_ind,gamm_gdpdef,del_gdpdef,...
    zstar,rstar,rkstar,wstar,wl_c,cstar,kstar,kbarstar,istar,ystar,sprd,zeta_spb,gammstar,vstar,nstar,...
    zeta_nRk,zeta_nR,zeta_nsigw,zeta_spsigw,zeta_nmue,zeta_spmue,zeta_nqk,zeta_nn] = getpara00_990(para)


nantpad = 20;
  
alp = para(1);
zeta_p = para(2);
iota_p = para(3);
del = .025;
ups = 1; %exp(para(4)/100);%maybe take out para(5)fix 0
Bigphi = para(5);
s2 = para(6);
h = para(7);
ppsi = para(8);
nu_l = para(9);
zeta_w = para(10);
iota_w = para(11);
law = 1.5;
laf = [];
bet  = 1/(1+para(12)/100);
psi1 = para(13);
psi2 = para(14);
psi3 = para(15);
pistar = (para(16)/100)+1;
sigmac = para(17);
rho = para(18);
epsp = 10;
epsw = 10;

%Financial Frictions Parameters
Fom = 1-(1-para(19))^(1/4);  %% F(omega) from annualized to quarterly default prob
sprd = (1+para(20)/100)^(1/4); %exp(para(21)/400);    %% st st spread from annual perc to quarterly number
zeta_spb = para(21);
gammstar = para(22);

npara = 22; 


%% exogenous processes - level

gam = para(npara+1)/100;
Lmean = para(npara+2);
gstar = .18;

npara = npara+2;

%% exogenous processes - autocorrelation

rho_g = para(npara+1);
rho_b = para(npara+2);
rho_mu = para(npara+3);
rho_z = para(npara+4);
rho_laf = para(npara+5);
rho_law = para(npara+6);
rho_rm = para(npara+7);

rho_sigw = para(npara+8);
rho_mue = para(npara+9);
rho_gamm = para(npara+10);

rho_pist = para(npara+11);
rho_lr = para(npara+12);

rho_zp = para(npara+13);
rho_tfp =para(npara+14);

rho_gdpdef = para(npara+15);
rho_pce    = para(npara+16);
npara = npara+16;


%% exogenous processes - standard deviation

sig_g = para(npara+1);
sig_b = para(npara+2);
sig_mu = para(npara+3);
sig_z = para(npara+4);
sig_laf = para(npara+5);
sig_law = para(npara+6);
sig_rm = para(npara+7);

sig_sigw = para(npara+8);
sig_mue = para(npara+9);
sig_gamm = para(npara+10);

sig_pist = para(npara+11);
sig_lr = para(npara+12);

sig_zp = para(npara+13);
sig_tfp= para(npara+14);

sig_gdpdef= para(npara+15);
sig_pce= para(npara+16);
npara = npara+16;

%% Standard deviations of the anticipated policy shocks
for i = 1:nantpad
  eval(strcat('sig_rm',num2str(i), '= para(npara +',num2str(i),');'));
  eval(strcat('sig_rm_ant(', num2str(i),') = sig_rm',num2str(i),';'));
end
npara = npara+nantpad;

eta_gz = para(npara+1);
eta_laf = para(npara+2);
eta_law = para(npara+3);

npara = npara+3;

modelalp_ind = para(npara+1);
gamm_gdpdef = para(npara+2);
del_gdpdef = para(npara+3);

npara = npara+3;


%% Parameters (implicit) -- from steady state

zstar = log(gam+1)+(alp/(1-alp))*log(ups); 

rstar = (1/bet)*exp(sigmac*zstar);

Rstarn = 100*(rstar*pistar-1);

rkstar = sprd*rstar*ups - (1-del);%NOTE: includes "sprd*"

wstar = (alp^(alp)*(1-alp)^(1-alp)*rkstar^(-alp)/Bigphi)^(1/(1-alp));

% Lstar = 1;
Lstar = (wstar/law/((1-gstar)*(alp/(1-alp)*wstar/rkstar)^alp/Bigphi-...
    (1-(1-del)/ups*exp(-zstar))*ups*exp(zstar)*alp/(1-alp)*wstar/rkstar)/(1-h*exp(-zstar)))^(1/(1+nu_l));

kstar = (alp/(1-alp))*wstar*Lstar/rkstar;

kbarstar = kstar*(gam+1)*ups^(1/(1-alp));

istar = kbarstar*( 1-((1-del)/((gam+1)*ups^(1/(1-alp)))) );

ystar = (kstar^alp)*(Lstar^(1-alp))/Bigphi;
if ystar <= 0

    disp([alp,  bet, kstar,Lstar])
    dm([ystar,Lstar,kstar,Bigphi])

end
cstar = (1-gstar)*ystar - istar;

wl_c = (wstar*Lstar/cstar)/law;


%FINANCIAL FRICTIONS ADDITIONS
% solve for sigmaomegastar and zomegastar
zwstar = norminv(Fom);
sigwstar = fzero(@(sigma)zetaspbfcn(zwstar,sigma,sprd)-zeta_spb,0.5);
% zetaspbfcn(zwstar,sigwstar,sprd)-zeta_spb % check solution

% evaluate omegabarstar
omegabarstar = omegafcn(zwstar,sigwstar);

% evaluate all BGG function elasticities
Gstar = Gfcn(zwstar,sigwstar);
Gammastar = Gammafcn(zwstar,sigwstar);
dGdomegastar = dGdomegafcn(zwstar,sigwstar);
d2Gdomega2star = d2Gdomega2fcn(zwstar,sigwstar);
dGammadomegastar = dGammadomegafcn(zwstar);
d2Gammadomega2star = d2Gammadomega2fcn(zwstar,sigwstar);
dGdsigmastar = dGdsigmafcn(zwstar,sigwstar);
d2Gdomegadsigmastar = d2Gdomegadsigmafcn(zwstar,sigwstar);
dGammadsigmastar = dGammadsigmafcn(zwstar,sigwstar);
d2Gammadomegadsigmastar = d2Gammadomegadsigmafcn(zwstar,sigwstar);

% evaluate mu, nk, and Rhostar
muestar = mufcn(zwstar,sigwstar,sprd);
nkstar = nkfcn(zwstar,sigwstar,sprd);
Rhostar = 1/nkstar-1;

% evaluate wekstar and vkstar
betbar=bet*exp((1-sigmac)*zstar);
wekstar = (1-gammstar/betbar)*nkstar...
    -gammstar/betbar*(sprd*(1-muestar*Gstar)-1);
vkstar = (nkstar-wekstar)/gammstar;

% evaluate nstar and vstar
nstar = nkstar*kbarstar;
vstar = vkstar*kbarstar;

% a couple of combinations
GammamuG = Gammastar-muestar*Gstar;
GammamuGprime = dGammadomegastar-muestar*dGdomegastar;

% elasticities wrt omegabar
zeta_bw = zetabomegafcn(zwstar,sigwstar,sprd);
zeta_zw = zetazomegafcn(zwstar,sigwstar,sprd);
zeta_bw_zw = zeta_bw/zeta_zw;

% elasticities wrt sigw
zeta_bsigw = sigwstar*(((1-muestar*dGdsigmastar/dGammadsigmastar)/...
    (1-muestar*dGdomegastar/dGammadomegastar)-1)*dGammadsigmastar*sprd+...
    muestar*nkstar*(dGdomegastar*d2Gammadomegadsigmastar-dGammadomegastar*d2Gdomegadsigmastar)/...
    GammamuGprime^2)/...
    ((1-Gammastar)*sprd+dGammadomegastar/GammamuGprime*(1-nkstar));
zeta_zsigw = sigwstar*(dGammadsigmastar-muestar*dGdsigmastar)/GammamuG;
zeta_spsigw = (zeta_bw_zw*zeta_zsigw-zeta_bsigw)/(1-zeta_bw_zw);

% elasticities wrt mue
zeta_bmue = -muestar*(nkstar*dGammadomegastar*dGdomegastar/GammamuGprime+dGammadomegastar*Gstar*sprd)/...
    ((1-Gammastar)*GammamuGprime*sprd+dGammadomegastar*(1-nkstar));
zeta_zmue = -muestar*Gstar/GammamuG;
zeta_spmue = (zeta_bw_zw*zeta_zmue-zeta_bmue)/(1-zeta_bw_zw);

% some ratios/elasticities
Rkstar = sprd*pistar*rstar; % (rkstar+1-delta)/ups*pistar;
zeta_Gw = dGdomegastar/Gstar*omegabarstar;
zeta_Gsigw = dGdsigmastar/Gstar*sigwstar;

% elasticities for the net worth evolution
zeta_nRk = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*(1-muestar*Gstar*(1-zeta_Gw/zeta_zw));
zeta_nR = gammstar/betbar*(1+Rhostar)*(1-nkstar+muestar*Gstar*sprd*zeta_Gw/zeta_zw);
zeta_nqk = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*(1-muestar*Gstar*(1+zeta_Gw/zeta_zw/Rhostar))...
    -gammstar/betbar*(1+Rhostar);
zeta_nn = gammstar/betbar+gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*zeta_Gw/zeta_zw/Rhostar;
zeta_nmue = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*(1-zeta_Gw*zeta_zmue/zeta_zw);
zeta_nsigw = gammstar*Rkstar/pistar/exp(zstar)*(1+Rhostar)*muestar*Gstar*(zeta_Gsigw-zeta_Gw/zeta_zw*zeta_zsigw);


end

function f=zetaspbfcn(z,sigma,sprd)
zetaratio = zetabomegafcn(z,sigma,sprd)/zetazomegafcn(z,sigma,sprd);
nk = nkfcn(z,sigma,sprd);
f = -zetaratio/(1-zetaratio)*nk/(1-nk);
end

function f=zetabomegafcn(z,sigma,sprd)
nk = nkfcn(z,sigma,sprd);
mustar = mufcn(z,sigma,sprd);
omegastar = omegafcn(z,sigma);
Gammastar = Gammafcn(z,sigma);
Gstar = Gfcn(z,sigma);
dGammadomegastar = dGammadomegafcn(z);
dGdomegastar = dGdomegafcn(z,sigma);
d2Gammadomega2star = d2Gammadomega2fcn(z,sigma);
d2Gdomega2star = d2Gdomega2fcn(z,sigma);
f = omegastar*mustar*nk*(d2Gammadomega2star*dGdomegastar-d2Gdomega2star*dGammadomegastar)/...
    (dGammadomegastar-mustar*dGdomegastar)^2/sprd/...
    (1-Gammastar+dGammadomegastar*(Gammastar-mustar*Gstar)/(dGammadomegastar-mustar*dGdomegastar));
end

function f=zetazomegafcn(z,sigma,sprd)
mustar = mufcn(z,sigma,sprd);
f = omegafcn(z,sigma)*(dGammadomegafcn(z)-mustar*dGdomegafcn(z,sigma))/...
    (Gammafcn(z,sigma)-mustar*Gfcn(z,sigma));
end

function f=nkfcn(z,sigma,sprd)
f = 1-(Gammafcn(z,sigma)-mufcn(z,sigma,sprd)*Gfcn(z,sigma))*sprd;
end

function f=mufcn(z,sigma,sprd)
f = (1-1/sprd)/(dGdomegafcn(z,sigma)/dGammadomegafcn(z)*(1-Gammafcn(z,sigma))+Gfcn(z,sigma));
end

function f=omegafcn(z,sigma)
f = exp(sigma*z-1/2*sigma^2);
end

function f=Gfcn(z,sigma)
f = normcdf(z-sigma);
end

function f=Gammafcn(z,sigma)
f = omegafcn(z,sigma)*(1-normcdf(z))+normcdf(z-sigma);
end

function f=dGdomegafcn(z,sigma)
f=normpdf(z)/sigma;
end

function f=d2Gdomega2fcn(z,sigma)
f = -z*normpdf(z)/omegafcn(z,sigma)/sigma^2;
end

function f=dGammadomegafcn(z)
f = 1-normcdf(z);
end

function f=d2Gammadomega2fcn(z,sigma)
f = -normpdf(z)/omegafcn(z,sigma)/sigma;
end

function f=dGdsigmafcn(z,sigma)
f = -z*normpdf(z-sigma)/sigma;
end

function f=d2Gdomegadsigmafcn(z,sigma)
f = -normpdf(z)*(1-z*(z-sigma))/sigma^2;
end

function f=dGammadsigmafcn(z,sigma)
f = -normpdf(z-sigma);
end

function f=d2Gammadomegadsigmafcn(z,sigma)
f = (z/sigma-1)*normpdf(z);
end
