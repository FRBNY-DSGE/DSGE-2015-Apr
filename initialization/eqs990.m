    %% Equilibrium conditions
euler     = 1;
inv       = 2;
capval    = 3;
spread    = 4;
nevol     = 5;
output    = 6;
caputl    = 7;
capsrv    = 8;
capev     = 9;
mkupp     = 10;
phlps     = 11;
caprnt    = 12;
msub      = 13;
wage      = 14;
mp        = 15;
res       = 16;
eq_g      = 17;
eq_b      = 18;
eq_mu     = 19;
eq_z      = 20;
eq_laf    = 21;
eq_law    = 22;
eq_rm     = 23;
eq_sigw   = 24;
eq_mue    = 25;
eq_gamm   = 26;
eq_laf1   = 27;
eq_law1   = 28;
eq_Ec     = 29;
eq_Eqk    = 30;
eq_Ei     = 31;
eq_Epi    = 32;
eq_EL     = 33;
eq_Erk    = 34;
eq_Ew     = 35;
eq_ERktil = 36;
euler_f   = 37;
inv_f     = 38;
capval_f  = 39;
output_f  = 40;
caputl_f  = 41;
capsrv_f  = 42;
capev_f   = 43;
mkupp_f   = 44;
caprnt_f  = 45;
msub_f    = 46;
res_f     = 47;
eq_Ec_f   = 48;
eq_Eqk_f  = 49;
eq_Ei_f   = 50;
eq_EL_f   = 51;
eq_Erk_f  = 52;
eq_ztil   = 53;
eq_pist   = 54;
pi1       = 55;
pi2       = 56;
pi_a      = 57;
Rt1       = 58;
eq_zp     = 59;
eq_Ez     = 60;

n_eqc     = 60;

if exist('nant','var')
  if nant > 0

    % These are the anticipated shocks. For each there is both an innovation
    % (for new anticipated shocks, calculated in period T only),
    % and a process, so that the shocks can be passed from period to
    % period.

    for i = 1:nant
      eval(strcat('eq_rml',num2str(i),'  = ',num2str(n_eqc+i),';'));
    end

    n_eqc=n_eqc+nant;

  end  
  
end
