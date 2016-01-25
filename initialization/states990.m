%% Endogenous variables 
y_t      = 1;
c_t      = 2;
i_t      = 3;
qk_t     = 4;
k_t      = 5;
kbar_t   = 6;
u_t      = 7;
rk_t     = 8;
Rktil_t  = 9;
n_t      = 10;
mc_t     = 11;
pi_t     = 12;
muw_t    = 13;%exo
w_t      = 14;
L_t      = 15;
R_t      = 16;
g_t      = 17;%exo
b_t      = 18;%exo
mu_t     = 19;%exo
z_t      = 20;%exo %%%z_t
laf_t    = 21;%exo
laf_t1   = 22;%exo
law_t    = 23;%exo
law_t1   = 24;%exo
rm_t     = 25;%exo
sigw_t   = 26;%exo - FF
mue_t    = 27;%exo - FF
gamm_t   = 28;%exo - FF
pist_t   = 29;%exo - pistart_t
E_c      = 30;%exp
E_qk     = 31;%exp
E_i      = 32;%exp
E_pi     = 33;%exp
E_L      = 34;%exp
E_rk     = 35;%exp
E_w      = 36;%exp
E_Rktil  = 37;%exp - FF
y_f_t    = 38;
c_f_t    = 39;
i_f_t    = 40;
qk_f_t   = 41;
k_f_t    = 42;
kbar_f_t = 43;
u_f_t    = 44;
rk_f_t   = 45;
w_f_t    = 46;
L_f_t    = 47;
r_f_t    = 48;
E_c_f    = 49;%exp
E_qk_f   = 50;%exp
E_i_f    = 51;%exp
E_L_f    = 52;%exp
E_rk_f   = 53;%exp4
ztil_t   = 54;
pi_t1    = 55;
pi_t2    = 56;
pi_a_t   = 57;
R_t1     = 58;
zp_t     = 59;
E_z      = 60;

nstates=60;


%/* shock indices */ EXOGENOUS
g_sh       = 1;
b_sh       = 2;
mu_sh      = 3;
z_sh       = 4;
laf_sh     = 5;
law_sh     = 6;
rm_sh      = 7;
sigw_sh    = 8;% - FF
mue_sh     = 9;% - FF
gamm_sh    = 10;% - FF
pist_sh    = 11;%
lr_sh      = 12; % Measurement error on the long rate

zp_sh      = 13;
tfp_sh     = 14;
gdpdef_sh  = 15;
pce_sh     = 16;

nex=16;

%/* expectation errors */
Ec_sh       = 1;
Eqk_sh      = 2;
Ei_sh       = 3;
Epi_sh      = 4;
EL_sh       = 5;
Erk_sh      = 6;
Ew_sh       = 7;
ERktil_sh   = 8;
Ec_f_sh     = 9;
Eqk_f_sh    = 10;
Ei_f_sh     = 11;
EL_f_sh     = 12;
Erk_f_sh    = 13;


n_exp=13;
n_exo=13;%12
n_end=nstates-n_exp-n_exo;
nend=n_exp;

if exist('nant','var')
  if nant > 0

    % These are the anticipated shocks. For each there is both an innovation
    % (for new anticipated shocks, calculated in period T only),
    % and a process, so that the shocks can be passed from period to
    % period.

    for i = 1:nant
      eval(strcat('rm_shl',num2str(i),' = ',num2str(nex + i),';'));
      eval(strcat('rm_tl',num2str(i),'  = ',num2str(nstates+i),';'));
    end

    n_exo = n_exo + nant;
    nex = nex + nant; 
    nstates=nstates+nant;

  end  

  
end

