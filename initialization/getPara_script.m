%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created 3/17/2014 by Matt Cocci

%% Summary/Motivation
% This script, depending on the mspec will assign values to the parameters
%   used in the mode
% This is done to allow for easy access to parameter values given a
%   parameter vector, rather than copyying the long list of output, as we
%   have below

%% Important Variables
% Important variables that must be set to run this:
%   1. para -- vector which will be assigned out to parameter names
%   2. mspec -- the mspec to use, which is crucial because different models
%       have different parameter names and numbers of parameters

%% Getting Parameter Values
% To access the modal values of the parameters
%   1. run your spec file
%   2. run forecast_mode_est_ant
%   3. set para = params
%   4. run this script

% Where this script is run:
%   1. dsgesolv.m
%   2. measure/meausurMSPEC.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if mspec==557
    [alp,zeta_p,iota_p,del,ups,Bigphi,s2,h,a2,nu_l,nu_m,zeta_w,iota_w,law,rstar,psi1,psi2,rho_r,pistar,...
        Fom,sprd,zeta_spb,gammstar,NEW_5,...
        gam,Lstar,chi,laf,gstar,Ladj,...
        rho_z,rho_phi,rho_chi,rho_laf,rho_mu,rho_b,rho_g,rho_sigw,rho_mue,rho_gamm,...
        sig_z,sig_phi,sig_chi,sig_laf,sig_mu,sig_b,sig_g,sig_r,sig_sigw,sig_mue,sig_gamm,...
        bet,zstar,phi,istokbarst,rkstar,cstar,ystar,istar,kstar,kbarstar,Rstarn,wstar,wadj,mstar,...
        zeta_nRk, zeta_nR, zeta_nqk, zeta_nn, zeta_nmue, zeta_spmue, zeta_nsigw, zeta_spsigw, ...
        vstar, nstar, Rkstar, sig_r_ant] = getpara00_557(para);
elseif mspec==990
    [alp,zeta_p,iota_p,del,ups,Bigphi,s2,h,ppsi,nu_l,zeta_w,iota_w,law,laf,bet,Rstarn,psi1,psi2,psi3,pistar,sigmac,rho,epsp,epsw...
        gam,Lmean,Lstar,gstar,rho_g,rho_b,rho_mu,rho_z,rho_laf,rho_law,rho_rm,rho_sigw,rho_mue,rho_gamm,rho_pist,rho_lr,rho_zp,rho_tfp,rho_gdpdef,rho_pce,...
        sig_g,sig_b,sig_mu,sig_z,sig_laf,sig_law,sig_rm,sig_sigw,sig_mue,sig_gamm,sig_pist,sig_rm_ant,sig_lr,sig_zp,sig_tfp,sig_gdpdef,sig_pce,...
        eta_gz,eta_laf,eta_law,modelalp_ind,gamm_gdpdef,del_gdpdef,...
        zstar,rstar,rkstar,wstar,wl_c,cstar,kstar,kbarstar,istar,ystar,sprd,zeta_spb,gammstar,vstar,nstar,...
        zeta_nRk,zeta_nR,zeta_nsigw,zeta_spsigw,zeta_nmue,zeta_spmue,zeta_nqk,zeta_nn] = getpara00_990(para);
end


