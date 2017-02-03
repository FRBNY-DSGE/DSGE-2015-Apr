%% Transformations:
%% format: [type, a, b, c]
%% Type 1:
%% x is [a,b] -> [-1,1] -> [-inf,inf] by (1/c)*c*z/sqrt(1-c*z^2)
%% Type 2:
%% x is [0,inf] -> [-inf,inf] by b + (1/c)*ln(para[i]-a);

function trspec = transp990(mspec)

	trspec = zeros(100,4);
  nantpad = 20;
	
	
	trspec(1,:) = [1	1E-5	.999	1]; 	%% alp;
	trspec(2,:) = [1  1E-5    0.999   1];   	%% zeta_p;
	trspec(3,:) = [1	1E-5	.999	1];	%% iota_p;
	trspec(4,:) = [2	1E-5	0	1];	%% ups;
	trspec(5,:) = [2	1.00	10.00	1];	%% Bigphi;
	trspec(6,:) = [0  -15.0   15.0    1];	%% s2;
	trspec(7,:) = [1	1E-5	.999	1]; 	%% h;
	trspec(8,:) = [1	1E-5	.999	1];	%% ppsi;
	trspec(9,:) = [2        1E-5     10     1];	%% nu_l;
	trspec(10,:) = [1       1E-5    0.999   1];   	%% zeta_w;
	trspec(11,:) = [1	1E-5	.999	1];	%% iota_w;
	trspec(12,:) = [2	1E-5	10	1];	%% bet;
	trspec(13,:) = [2	1E-5	10.00	1];	%% psi1;
	trspec(14,:) = [0	-0.5	0.5	1];	%% psi2;
  trspec(15,:) = [0	-0.5	0.5	1];	%% psi3;
	trspec(16,:) = [2	1E-5	10	1];	%% pistar;
  trspec(17,:) = [2	1E-5	10	1];	%% sigmac;
  trspec(18,:) = [1	1E-5	.999	1];	%% rho;

  trspec(19,:) = [1	1E-5	.99	1];     %% F(omega)
	trspec(20,:) = [2	1E-5	0	1];     %% st st spread
  trspec(21,:) = [1	1E-5	.99	1];     %% zeta_sp
  trspec(22,:) = [1	1E-5	.99	1];     %% gammst

                
	npara = 22;
	
	%% exogenous processes - level
	trspec(npara+1,:) = [0	-5.0	5.0	1];	%% gam;
  	trspec(npara+2,:) = [0   -1000    1000      1]; %% Lmean;
	
	npara = npara+2;	
	
	%% exogenous processes - autocorrelation
	trspec(npara+1,:) = [1	1E-5	.999	1];	%% rho_g;
	trspec(npara+2,:) = [1	1E-5	.999	1];	%% rho_b;
	trspec(npara+3,:) = [1	1E-5	.999	1];	%% rho_mu;
	trspec(npara+4,:) = [1	1E-5	.999	1];	%% rho_z;
	trspec(npara+5,:) = [1	1E-5	.999	1];	%% rho_laf;
	trspec(npara+6,:) = [1	1E-5	.999	1];	%% rho_law;
	trspec(npara+7,:) = [1	1E-5	.999	1];	%% rho_rm;
  trspec(npara+8,:) = [1	1E-5	.999	1];	%% rho_sigw;
  trspec(npara+9,:) = [1	1E-5	.99	1];	%% rho_mue;
  trspec(npara+10,:) = [1	1E-5	.99	1];	%% rho_gamm;
	trspec(npara+11,:) = [1	1E-5	.999	1];	%% rho_pist;
	trspec(npara+12,:) = [1	1E-5	.999	1];	%% rho_lr;
        
	trspec(npara+13,:) = [1	1E-5	.999	1];	%% rho_zp;
  trspec(npara+14,:) = [1	1E-5	.999	1];	%% rho_tfp;
  trspec(npara+15,:) = [1	1E-5	.999	1];	%% rho_gdpdef;
  trspec(npara+16,:) = [1	1E-5	.999	1];	%% rho_pce;

	npara = npara+16;	
    
	%% exogenous processes - standard deviation
	trspec(npara+1,:) = [2	1E-8	5	1];	%% sig_g;
	trspec(npara+2,:) = [2	1E-8	5	1];	%% sig_b;
	trspec(npara+3,:) = [2	1E-8	5	1];	%% sig_mu;
	trspec(npara+4,:) = [2	1E-8	5	1];	%% sig_z;
	trspec(npara+5,:) = [2	1E-8	5	1];	%% sig_laf;
	trspec(npara+6,:) = [2	1E-8	5	1];	%% sig_law;
	trspec(npara+7,:) = [2	1E-8	5	1];	%% sig_rm;
  trspec(npara+8,:) = [2	1E-5	0	1];	%% sig_sigw;
  trspec(npara+9,:) = [2	1E-5	0	1];	%% sig_mue;
  trspec(npara+10,:) = [2	1E-5	0	1];	%% sig_gamm;
	trspec(npara+11,:) = [2	1E-8	5	1];	%% sig_pist;        
	trspec(npara+12,:) = [2	1E-8	5	1];	%% sig_lr;        

	trspec(npara+13,:) = [2	1E-8	5	1];	%% sig_zp;
  trspec(npara+14,:) = [2	1E-8	5	1];	%% sig_tfp;

  trspec(npara+15,:) = [2	1E-8	5	1];	%% sig_gdpdef;
  trspec(npara+16,:) = [2	1E-8	5	1];	%% sig_pce;
	npara = npara+16;
        
  %% Standard Deviation of the Anticipated Shocks
  for i = 1:nantpad
      eval(strcat('trspec(npara +',num2str(i),',:) = [2 1E-5  0  1];'));
  end
  npara = npara+nantpad;

  trspec(npara+1,:) = [1   1E-5  0.999  1];	%% eta_gz;
  trspec(npara+2,:) = [1   1E-5  0.999  1];	%% eta_laf;
  trspec(npara+3,:) = [1   1E-5  0.999  1];	%% eta_law;       
	npara = npara+3;

  trspec(npara+1,:) = [0   0     0      0];	%% modelalp_ind;
  trspec(npara+2,:) = [0  -10   -10     1];	%% gamm_gdpdef;
  trspec(npara+3,:) = [0  -10   -10     1];	%% del_gdpdef;

	npara = npara+3;


	trspec = trspec(1:npara,:);



