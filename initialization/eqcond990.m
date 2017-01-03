%%*****************************************************
%**      1. Consumption Euler Equation
%******************************************************/
%%* sticky prices and wages */
G0(euler,c_t) =  1;
G0(euler,R_t) = (1-h*exp(-zstar))/(sigmac*(1+h*exp(-zstar)));
G0(euler,b_t) = -1;
G0(euler,E_pi) = -(1-h*exp(-zstar))/(sigmac*(1+h*exp(-zstar)));
G0(euler,z_t) = (h*exp(-zstar))/(1+h*exp(-zstar));
G0(euler,E_c) = -1/(1+h*exp(-zstar));
%G0(euler,ztil_t) = -( 1/(1-alp) )*(rho_z-1)/(1+h*exp(-zstar));
G0(euler,E_z) = -1/(1+h*exp(-zstar));
G0(euler,L_t) = -(sigmac - 1)*wl_c/(sigmac*(1 + h*exp(-zstar)));
G0(euler,E_L) = (sigmac - 1)*wl_c/(sigmac*(1 + h*exp(-zstar)));
G1(euler,c_t) = (h*exp(-zstar))/(1+h*exp(-zstar));

%%* flexible prices and wages **/
G0(euler_f,c_f_t) =  1;
G0(euler_f,r_f_t) = (1-h*exp(-zstar))/(sigmac*(1+h*exp(-zstar)));
G0(euler_f,b_t) = -1;
G0(euler_f,z_t) =   (h*exp(-zstar))/(1+h*exp(-zstar));
G0(euler_f,E_c_f) = -1/(1+h*exp(-zstar));
%G0(euler_f,ztil_t) = -( 1/(1-alp) )*(rho_z-1)/(1+h*exp(-zstar));
G0(euler_f,E_z) = -1/(1+h*exp(-zstar));
G0(euler_f,L_f_t) = -(sigmac - 1)*wl_c/(sigmac*(1 + h*exp(-zstar)));
G0(euler_f,E_L_f) = (sigmac - 1)*wl_c/(sigmac*(1 + h*exp(-zstar)));
G1(euler_f,c_f_t) = (h*exp(-zstar))/(1+h*exp(-zstar));

%%****************************************************
%**      2. Investment Euler Equation
%*****************************************************/
%%* sticks prices and wages **/
G0(inv,qk_t) = -1/(s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar)));
G0(inv,i_t) = 1;
G0(inv,z_t) = 1/(1+bet*exp((1-sigmac)*zstar));
G1(inv,i_t) = 1/(1+bet*exp((1-sigmac)*zstar));
G0(inv,E_i) = -bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
%G0(inv,ztil_t) = -( 1/(1-alp) )*(rho_z-1)*bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
G0(inv,E_z) = -bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
G0(inv,mu_t) = -1;

%%* flexible prices and wages **/
G0(inv_f,qk_f_t) = -1/(s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar)));
G0(inv_f,i_f_t) = 1;
G0(inv_f,z_t) = 1/(1+bet*exp((1-sigmac)*zstar));
G1(inv_f,i_f_t) = 1/(1+bet*exp((1-sigmac)*zstar));
G0(inv_f,E_i_f) = -bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
%G0(inv_f,ztil_t) = -( 1/(1-alp) )*(rho_z-1)*bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
G0(inv_f,E_z) = -bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar));
G0(inv_f,mu_t) = -1;

%%****************************************************
%**      3. FINANCIAL FRICTION BLOCK
%*****************************************************/

%% return to capital
%%* sticky prices and wages **/
G0(capval,Rktil_t) = 1;
G0(capval,pi_t) = -1;
G0(capval,rk_t) = -rkstar/(rkstar+1-del);
G0(capval,qk_t) = -(1-del)/(rkstar+1-del);
G1(capval,qk_t) = -1;

%% spreads
%%* sticky prices and wages **/
G0(spread,E_Rktil) = 1;
G0(spread,R_t) = -1;
G0(spread,b_t) = (sigmac*(1+h*exp(-zstar)))/(1-h*exp(-zstar));
G0(spread,qk_t) = -zeta_spb;
G0(spread,kbar_t) = -zeta_spb;
G0(spread,n_t) = zeta_spb;
G0(spread,sigw_t) = -1;
G0(spread,mue_t) = -1;

%% n evol
%%* sticky prices and wages **/
G0(nevol,n_t) = 1;
G0(nevol,gamm_t) = -1;
G0(nevol,z_t) = gammstar*vstar/nstar;
G0(nevol,Rktil_t) = -zeta_nRk;
G0(nevol,pi_t) = (zeta_nRk - zeta_nR);
G1(nevol,sigw_t) = -zeta_nsigw/zeta_spsigw;
G1(nevol,mue_t) = -zeta_nmue/zeta_spmue;
G1(nevol,qk_t) = zeta_nqk;
G1(nevol,kbar_t) = zeta_nqk;
G1(nevol,n_t) = zeta_nn;
G1(nevol,R_t) = -zeta_nR;%BUG FIXED -- TERM ADDED ON 9/27/11
G1(nevol,b_t) = zeta_nR*((sigmac*(1+h*exp(-zstar)))/(1-h*exp(-zstar)));% BUG FIXED -- TERM NORMALIZED ON 01/03/17

%%* flexible prices and wages - ASSUME NO FINANCIAL FRICTIONS
G0(capval_f,E_rk_f) = -rkstar/(rkstar+1-del);
G0(capval_f,E_qk_f) = -(1-del)/(rkstar+1-del);
G0(capval_f,qk_f_t) = 1;
G0(capval_f,r_f_t) = 1;
G0(capval_f,b_t) = -(sigmac*(1+h*exp(-zstar)))/(1-h*exp(-zstar));


%%***************************************************
%**      4. Aggregate Production Function
%****************************************************/
%%* sticky prices and wages **/
G0(output,y_t ) =  1;
G0(output,k_t) = -Bigphi*alp;
G0(output,L_t) = -Bigphi*(1-alp);
%G0(output,ztil_t) = - ( Bigphi-1 ) / (1-alp);  % See supersticky section on adding long run changes to productivity

%%* flexible prices and wages **/
G0(output_f,y_f_t ) =  1;
G0(output_f,k_f_t) = -Bigphi*alp;
G0(output_f,L_f_t) = -Bigphi*(1-alp);
%G0(output_f,ztil_t) = - ( Bigphi-1 ) / (1-alp); % See supersticky section on adding long run changes to productivity

%%**************************************************
%**      5. Capital Utilization
%***************************************************/
%%* sticky prices and wages **/
G0(caputl,k_t) =  1;
G1(caputl,kbar_t ) =  1;
G0(caputl,z_t ) = 1;
G0(caputl,u_t) = -1;

%%* flexible prices and wages **/
G0(caputl_f,k_f_t) =  1;
G1(caputl_f,kbar_f_t ) =  1;
G0(caputl_f,z_t ) = 1;
G0(caputl_f,u_f_t) = -1;

%%*************************************************
%**      6. Rental Rate of Capital
%**************************************************/
%%* sticky prices and wages **/
G0(capsrv,u_t ) = 1;
G0(capsrv,rk_t) = -(1-ppsi)/ppsi;

%% flexible prices and wages **/
G0(capsrv_f,u_f_t ) = 1;
G0(capsrv_f,rk_f_t) = -(1-ppsi)/ppsi;

%%*************************************************
%**      7. Evolution of Capital
%**************************************************/
%% sticky prices and wages **/
G0(capev,kbar_t) =  1;
G1(capev,kbar_t) =  1-istar/kbarstar;
G0(capev,z_t) =  1-istar/kbarstar;
G0(capev,i_t) = -istar/kbarstar;
G0(capev,mu_t) = -istar*s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar))/kbarstar;

%% flexible prices and wages **/
G0(capev_f,kbar_f_t) =  1;
G1(capev_f,kbar_f_t) =  1-istar/kbarstar;
G0(capev_f,z_t) =  1-istar/kbarstar;
G0(capev_f,i_f_t) = -istar/kbarstar;
G0(capev_f,mu_t) = -istar*s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar))/kbarstar;


%%***********************************************
%**      8. Price Markup
%************************************************/
%%* sticky prices and wages **/
G0(mkupp,mc_t) = 1;
G0(mkupp,w_t) =  -1;
G0(mkupp,L_t) =  -alp;
G0(mkupp,k_t) =  alp;

%%* flexible prices and wages **/
G0(mkupp_f,w_f_t) =  1;
G0(mkupp_f,L_f_t) =  alp;
G0(mkupp_f,k_f_t) =  -alp;

%%***********************************************
%**      9. Phillips Curve
%************************************************/
%%* sticky prices and wages **/
G0(phlps,pi_t) =  1;
G0(phlps,mc_t) =  -((1-zeta_p*bet*exp((1-sigmac)*zstar))*(1-zeta_p))/(zeta_p*((Bigphi-1)*epsp+1))*1/(1+iota_p*bet*exp((1-sigmac)*zstar));
G1(phlps,pi_t) = iota_p*1/(1+iota_p*bet*exp((1-sigmac)*zstar));
G0(phlps,E_pi) = -bet*exp((1-sigmac)*zstar)*1/(1+iota_p*bet*exp((1-sigmac)*zstar));
% Comment out for counterfactual with no price mark up shock
G0(phlps,laf_t) = -(1+iota_p*bet*exp((1-sigmac)*zstar))*1/(1+iota_p*bet*exp((1-sigmac)*zstar));

%%* flexible prices and wages **/
%%* not necesary **/

%%**********************************************
%**     10. Rental Rate of Capital
%***********************************************/
%%* sticky prices and wages **/
G0(caprnt,rk_t) =  1;
G0(caprnt,k_t) =  1;
G0(caprnt,L_t) =  -1;
G0(caprnt,w_t) =  -1;

%%* flexible prices and wages **/
G0(caprnt_f,rk_f_t) =  1;
G0(caprnt_f,k_f_t) =  1;
G0(caprnt_f,L_f_t) =  -1;
G0(caprnt_f,w_f_t) =  -1;


%%*********************************************
%**     11. Marginal Substitution
%***********************************************/
%%* sticky prices and wages **/
G0(msub,muw_t) =  1;
G0(msub,L_t) =  nu_l;
G0(msub,c_t) = 1/( 1-h*exp(-zstar) );
G1(msub,c_t) = h*exp(-zstar)/( 1-h*exp(-zstar) );
G0(msub,z_t) = h*exp(-zstar) /( 1-h*exp(-zstar) );
G0(msub,w_t) = -1;

%%* flexible prices and wages **/
G0(msub_f,w_f_t) =  -1;
G0(msub_f,L_f_t) =  nu_l;
G0(msub_f,c_f_t) = 1/( 1-h*exp(-zstar) );
G1(msub_f,c_f_t) = h*exp(-zstar)/( 1-h*exp(-zstar) );
G0(msub_f,z_t) = h*exp(-zstar)/( 1-h*exp(-zstar) );

%%********************************************
%**     12. Evolution of Wages
%**********************************************/
%%* sticky prices and wages **/
G0(wage,w_t)    = 1;
G0(wage,muw_t)  = (1-zeta_w*bet*exp((1-sigmac)*zstar))*(1-zeta_w)/(zeta_w*((law-1)*epsw+1))*1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,pi_t)   = (1+iota_w*bet*exp((1-sigmac)*zstar))*1/(1+bet*exp((1-sigmac)*zstar));
G1(wage,w_t)    = 1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,z_t)    = 1/(1+bet*exp((1-sigmac)*zstar));
G1(wage,pi_t)   = iota_w*1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,E_w)    = -bet*exp((1-sigmac)*zstar)*1/(1+bet*exp((1-sigmac)*zstar));
%G0(wage,ztil_t) = -( 1/(1-alp) )*(rho_z-1)*bet*exp((1-sigmac)*zstar)*1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,E_z) = -bet*exp((1-sigmac)*zstar)*1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,E_pi)   = -bet*exp((1-sigmac)*zstar)*1/(1+bet*exp((1-sigmac)*zstar));
G0(wage,law_t)  = -1;

%%* flexible prices and wages **/
%%* not necessary **/

%%*********************************************
%**    13. Monetary Policy Rule
%***********************************************/
%%* sticky prices and wages **/
G0(mp,R_t) = 1;
G1(mp,R_t) = rho;
G0(mp,pi_t) = -(1-rho)*psi1;
G0(mp,pist_t) = (1-rho)*psi1;
G0(mp,y_t) = -(1-rho)*psi2 - psi3;
G0(mp,y_f_t) = (1-rho)*psi2 + psi3;
%G0(mp,y_t) = - psi3;
%G0(mp,y_f_t) = psi3;
G1(mp,y_t) =  -psi3;
G1(mp,y_f_t) =  psi3;
G0(mp,rm_t )   = -1;

%%* flexible prices and wages **/
%%* not necessary **/

%%********************************************
%**   14. Resource Constraint
%***********************************************/
%%* sticky prices and wages **/
G0(res,y_t) = 1;
G0(res,g_t) = -gstar;
%if rho_z<1;
%  G0(res,ztil_t) = gstar*( 1/(1-alp) ); % See supersticky section on adding long run changes to productivity
%end;
G0(res,c_t) = -cstar/ystar;
G0(res,i_t) = -istar/ystar;
G0(res,u_t) = -rkstar*kstar/ystar;

%%* flexible prices and wages **/
G0(res_f,y_f_t) = 1;
G0(res_f,g_t) = -gstar;
%if rho_z<1;
%  G0(res_f,ztil_t) = gstar*( 1/(1-alp) );% See supersticky section on adding long run changes to productivity
%end;
G0(res_f,c_f_t) = -cstar/ystar;
G0(res_f,i_f_t) = -istar/ystar;
G0(res_f,u_f_t) = -rkstar*kstar/ystar;

%%*********************************************
%**   Exogenous Processes
%***********************************************/

%%* neutral technology **/
G0(eq_z,z_t) = 1;
G1(eq_z,ztil_t) = ( 1/(1-alp) )*(rho_z-1);
G0(eq_z,zp_t) = -1;
 PSI(eq_z,z_sh) = ( 1/(1-alp) );

G0(eq_ztil,ztil_t) = 1;
G1(eq_ztil,ztil_t) = rho_z;
 PSI(eq_ztil,z_sh) = 1;

% Extra term for long run changes to productivity; added by MDC
G0(eq_zp, zp_t) = 1;
G1(eq_zp, zp_t) = rho_zp;
PSI(eq_zp, zp_sh) = 1;

%%* government spending **/
G0(eq_g,g_t)  =  1;
G1(eq_g,g_t)  =  rho_g;
 PSI(eq_g,g_sh)  =  1;
 PSI(eq_g,z_sh)  =  eta_gz;

%%* asset shock **/
G0(eq_b,b_t)  =  1;
G1(eq_b,b_t)  =  rho_b;
 PSI(eq_b,b_sh)  =  1;

%%* investment specific technology **/
G0(eq_mu,mu_t) = 1;
G1(eq_mu,mu_t) = rho_mu;
 PSI(eq_mu,mu_sh) = 1;

%/** price mark-up shock **/
G0(eq_laf,laf_t)  =  1;
G1(eq_laf,laf_t)  =  rho_laf;
G1(eq_laf,laf_t1) =  -eta_laf;
 PSI(eq_laf,laf_sh)  =  1;

G0(eq_laf1,laf_t1) = 1;
 PSI(eq_laf1,laf_sh ) = 1;

%/** wage mark-up shock **/
G0(eq_law,law_t)  =  1;
G1(eq_law,law_t)  =  rho_law;
G1(eq_law,law_t1) =  -eta_law;
 PSI(eq_law,law_sh)  =  1;

G0(eq_law1,law_t1) = 1;
 PSI(eq_law1,law_sh ) = 1;

 %/** monetary policy shock **/
G0(eq_rm,rm_t) = 1;
G1(eq_rm,rm_t) = rho_rm;
 PSI(eq_rm,rm_sh) = 1;


%FINANCIAL FRICTIONS:

%% sigw shock
G0(eq_sigw,sigw_t) = 1;
G1(eq_sigw,sigw_t) = rho_sigw;
PSI(eq_sigw,sigw_sh) = 1;


%% mue shock
G0(eq_mue,mue_t) = 1;
G1(eq_mue,mue_t) = rho_mue;
PSI(eq_mue,mue_sh) = 1;


%% gamm shock
G0(eq_gamm,gamm_t) = 1;
G1(eq_gamm,gamm_t) = rho_gamm;
PSI(eq_gamm,gamm_sh) = 1;


%% long term inflation expectations
 %/** pistar **/
G0(eq_pist,pist_t) = 1;
G1(eq_pist,pist_t) = rho_pist;
 PSI(eq_pist,pist_sh) = 1;


%/** ANTICIPATED SHOCKS **/
if exist('nant','var')

  if nant > 0

    % This section adds the anticipated shocks. There is one state for all the
    % anticipated shocks that will hit in a given period (i.e. rm_tl2 holds
    % those that will hit in two periods), and the equations are set up so that
    % rm_tl2 last period will feed into rm_tl1 this period (and so on for other
    % numbers), and last period's rm_tl1 will feed into the rm_t process (and
    % affect the Taylor Rule this period).

    %note: belongs to eq_rm equation above^^
        G1(eq_rm,rm_tl1) = 1;
    %___________________________________

      G0(eq_rml1,rm_tl1) = 1;
      PSI(eq_rml1,rm_shl1) = 1;

      if nant > 1
        for i = 2:nant
          eval(strcat('G1(eq_rml',num2str(i-1),',rm_tl',num2str(i),') = 1;'));
          %___________________________________

          eval(strcat('G0(eq_rml',num2str(i),',rm_tl',num2str(i),') = 1;'));
          eval(strcat('PSI(eq_rml',num2str(i),',rm_shl',num2str(i),') = 1;'));
        end
      end
    end
end

%******************************************
%**    Rational Expectations Errors
%*******************************************/
%* E[c) **/
%* sticky prices and wages **/
G0(eq_Ec,c_t ) = 1;
G1(eq_Ec,E_c) = 1;
 PIE(eq_Ec,Ec_sh ) = 1;

%* flexible prices and wages **/
G0(eq_Ec_f,c_f_t ) = 1;
G1(eq_Ec_f,E_c_f) = 1;
 PIE(eq_Ec_f,Ec_f_sh ) = 1;

%* E(q) **/
%* sticky prices and wages **/
G0(eq_Eqk,qk_t) = 1;
G1(eq_Eqk,E_qk) = 1;
 PIE(eq_Eqk,Eqk_sh ) = 1;

%* flexible prices and wages **/
G0(eq_Eqk_f,qk_f_t ) = 1;
G1(eq_Eqk_f,E_qk_f) = 1;
 PIE(eq_Eqk_f,Eqk_f_sh) = 1;

%* E(i) **/
%* sticky prices and wages **/
G0(eq_Ei,i_t ) = 1;
G1(eq_Ei,E_i) = 1;
 PIE(eq_Ei,Ei_sh ) = 1;

%* flexible prices and wages **/
G0(eq_Ei_f,i_f_t ) = 1;
G1(eq_Ei_f,E_i_f) = 1;
 PIE(eq_Ei_f,Ei_f_sh ) = 1;

%* E(pi) **/
%* sticky prices and wages **/
G0(eq_Epi,pi_t ) = 1;
G1(eq_Epi,E_pi) = 1;
 PIE(eq_Epi,Epi_sh ) = 1;

%* E(l) **/
%* sticky prices and wages **/
G0(eq_EL,L_t ) = 1;
G1(eq_EL,E_L) = 1;
 PIE(eq_EL,EL_sh ) = 1;

%* flexible prices and wages **/
G0(eq_EL_f,L_f_t ) = 1;
G1(eq_EL_f,E_L_f) = 1;
 PIE(eq_EL_f,EL_f_sh ) = 1;

%* E(rk) **/
%* sticky prices and wages **/
G0(eq_Erk,rk_t ) = 1;
G1(eq_Erk,E_rk) = 1;
 PIE(eq_Erk,Erk_sh ) = 1;

%* flexible prices and wages **/
G0(eq_Erk_f,rk_f_t ) = 1;
G1(eq_Erk_f,E_rk_f) = 1;
 PIE(eq_Erk_f,Erk_f_sh ) = 1;

%* E(w) **/
%* sticky prices and wages **/
G0(eq_Ew,w_t ) = 1;
G1(eq_Ew,E_w) = 1;
PIE(eq_Ew,Ew_sh ) = 1;

%% E_Rktil
%* sticky prices and wages **/
G0(eq_ERktil,Rktil_t) = 1;
G1(eq_ERktil,E_Rktil) = 1;
PIE(eq_ERktil,ERktil_sh) = 1;


%******************************************
%% EXTRA STATES
% These aren't strictly necessary, but they
% track lags or simplify the equations
%*******************************************/


%% pi_t1
G0(pi1,pi_t1) = 1;
G1(pi1,pi_t) = 1;

%% pi_t2
G0(pi2,pi_t2) = 1;
G1(pi2,pi_t1) = 1;

%% pi_a
G0(pi_a, pi_a_t)=  1;
G0(pi_a, pi_t)  = -1;
G0(pi_a, pi_t1) = -1;
G0(pi_a, pi_t2) = -1;
G1(pi_a, pi_t2) =  1;

%% Rt1
G0(Rt1, R_t1)= 1;
G1(Rt1, R_t) = 1;

%% E_z
G0(eq_Ez, E_z)    = 1;
G0(eq_Ez, ztil_t) = -(1/(1-alp))*(rho_z-1);
G0(eq_Ez, zp_t)   = -rho_zp;
