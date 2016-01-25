%% gibb_est_ant.m
% This program produces and saves draws from the posterior distribution of
% the parameters.

%% Step 1: Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Basic Initialization (Reads in data and model specifications)
initializePrograms;
marglh = 'marglh';

%%  Configure the Metropolis-Hastings Algorithm

% Set jump size parameters for sampling
cc0 = 0.01;
cc = .09;
date_q = (1:1:qahead)';


%% Initialize Mode-Finding Algorithm

% Here, we detect whether we have already gone through gibb_est_ant.m. If
% not, we start the mode-finding algorithm from the pre-specified 'para' vector. 
% Otherwise, we will begin the optimization from where we left off in the last
% round of iterations.

if(~exist('xh', 'var'));

    % specify starting mode
    infile0 = [spath, 'mode_in'];
    fid0 = fopen(infile0,'r');
    params = fread(fid0,[npara,1],'single');
    disp('Previous mode found, reading in...')
else
    disp('gibb_est_ant was called recursively. Will use parameter vector from last call to find mode');
end

% Inputs to minimization algorithm (csminwel).
x0 = invtrans(params,trspec);
H0 = eye(npara)*1E-4;
nit = 1000;
crit= 1E-10;
MIN = 1;

%% Step 2: Find mode of posterior distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if reoptimize
  disp('Re-optimizing')
  % We seek the posterior distribution of the parameters.
  % We first find the mode of the distribution (i.e., the maximum of its pdf)
  % so that once we begin sampling, we can begin from this mode.

  [fh,xh,g,H,itct,fcount,retcode] = csminwel('objfcndsge_2part',x0,H0,[],crit,nit,...
      YY,nobs,nlags,nvar, mspec,npara,trspec,pmean,pstdd,pshape,para_mask,...
      para_fix,marglh,coint,cointadd,cointall,MIN, nant, antlags);

  % Transform modal parameters so they are no longer bounded (i.e., allowed 
  % to lie anywhere on the real line).

  xh = real(xh);
  params = trans(xh,trspec);

  % If we choose to fix any of the parameters (i.e., not estimate them), here
  % we reset those parameters to those fixed values.

  params = params.*(1-para_mask)+para_fix.*para_mask;

  % Once we get a mode, we this save the parameter vector in the binary file
  % 'outfile0'.

  outfile0 = [spath,'/mode_out'];
  fid0 = fopen(outfile0,'w');
  fwrite(fid0,params','single');
  fclose(fid0);

  % If the algorithm stops only because we have exceeded the maximum number of
  % iterations,continue improving guess of modal parameters by recursively
  % calling gibb.
  if itct >= nit
      clear infile0;
      gibb_est_ant;
      return;
  end
end


%% Step 3: Compute the inverse Hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Once we find and save the posterior mode, we calculate the Hessian
% at the mode. The hessian is used to calculate the variance of the proposal 
% distribution,which is used to draw a new parameter in each iteration of 
% the algorithm.
%% Step 3a: Calculate Hessian corresponding to the posterior mode
if CH == 1 % if choose to calculate Hessian
    disp('Re-computing Hessian')
    MIN=0;
    %calculate Hessian
    hessian = hessizero('objfcndsge_2part',[params,para_mask],1,...
        YY,nobs,nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,para_mask,...
        para_fix,marglh,coint,cointadd,cointall,MIN,nant,antlags);
    % Save computed Hessian in output file

    outfile_hes = [spath,'/hessian'];    
    fid_hes = fopen(outfile_hes,'w');
    fwrite(fid_hes,hessian,'single');
    fclose(fid_hes);
else % if don't calculate Hessian, read in pre-calculated Hessian matrix
    disp('Using pre-calculated Hessian')
    infile_hes = [spath,'/hessian']; 
    fid = fopen(infile_hes,'r');
    hessian=fread(fid,[npara npara],'single');
end;

if any(diag(hessian)<0)
    error('negative in diagonal of hessian');
end

%% Step 3b: Calculate inverse of Hessian (by Singular Value Decomposition)
% We use singular value decomposition to compute the inverse of the
% Hessian, as the Hessian calculated above is singular
rankhhm = npara;
[u,s,v] = svd(hessian);
sigpropinv  =  hessian;
md = min(size(s));
bigev = find(diag(s(1:md,1:md))>1e-6);
sigpropdim = length(bigev);

sigproplndet = 0;

for i = 1:npara
    if i > sigpropdim
        s(i,i) = 0;
    else
        s(i,i)    = 1/s(i,i);
        sigproplndet = sigproplndet+log(s(i,i));
    end
end

invhhm  = u*s*u';
sigscale = u*sqrt(s);

if length(bigev) ~= (npara-length(find(para_mask)))
    disp('problem -- shutting down dimensions')
    pause
end


%% Step 3: Initialize algorithm by drawing para_old from a normal distribution centered on the posterior mode (propdens).

% para_old is used to solve the model to check for indeterminancy, and it is
% passed through objfcnmhdsge.m to calculate the posterior value.

% objfcnmhdsge.m is similar to objfcndsge.m, except it also checks that parameters
% are within the bounds specified in mspec_parameters<mspec>.m. Until the parameters
% are within bounds, or until the posterior value is sufficiently large, gibb.m
% keeps drawing a new para_old (if you see multiple lines of "Initializing..."
% in the command window, this is the problem).

[TTT,RRR,CCC,valid] = dsgesolv(mspec,params,nant);
retcode = valid;
[lnpost0,lnpy0,zend0,ZZ0,DD0,QQ0] = feval('objfcnmhdsge_2part',params,bounds,YY,YY0,nobs,...
                                    nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,...
                                    TTT,RRR,CCC,valid,para_mask,coint,cointadd,cointall,YYcoint0,nant,antlags,nshocks);

jc = 1;
valid0 = 0;
while ~valid0
    jc = jc+1;
    
    para_old   = params + cc0*(sigscale*randn(npara,1));
    para_old = para_old.*(1-para_mask)+para_fix.*para_mask;
    
    [TTT_old,RRR_old,CCC_old,valid_old] = dsgesolv(mspec,para_old,nant);
    nstate = length(TTT_old);
    nshocks = size(RRR_old,2);
    
    retcode = valid_old;
    [post_old,like_old,zend_old,ZZ_old,DD_old,QQ_old] = feval('objfcnmhdsge_2part',para_old,bounds,YY,YY0,nobs,...
                                                        nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,...
                                                        TTT_old,RRR_old,CCC_old,valid_old,para_mask,coint,cointadd,cointall,YYcoint0,nant,antlags,nshocks-nant);
                                                    
    propdens  = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - 0.5*sigpropdim*log(cc0^2) ...
                -0.5*(para_old - params)'*sigpropinv*(para_old - params)/cc0^2;
    
    if post_old > -10000000;
        valid0 = 1;
    end
    
    record(jc) = post_old;
    disp('Initializing Metropolis-Hastings Algorithm')
    
end

tic
fprintf(1,'\n Time %3.2f ',ti(I));
fprintf(1,'\n Peak = %2.3f',lnpost0);

%% Open files for saving
% Once you open files for writing, the original files (if they existed)
% will be overwritten. 

% Parameters draws
outfile1 = [spath,'/params'];
fid1 = fopen(outfile1,'w');

% Transition Matrices
outfile2 = [spath,'/post'];
fid2 = fopen(outfile2,'w');

outfile3 = [spath,'/TTT'];
fid3 = fopen(outfile3,'w');

outfile4 = [spath,'/RRR'];
fid4 = fopen(outfile4,'w');

outfile5 = [spath,'/zend'];
fid5 = fopen(outfile5,'w');

% Covariance Matrices
outfile6 = [spath,'/cov'];
fid6 = fopen(outfile6,'w');

% Initialize some variables (for calculating rejection rate)
Tim = 0;
eT = 0;
reje = 0;

%% Step 4: For nsim*ntimes iterations within each block, generate a new parameter draw.
%%         Decide to accept or reject, and save every ntimes_th draw that is accepted.

for iblock = 1:nblocks
    
    if iblock >1; fprintf(1,' block: %2.0f %3.2f ;',[iblock,toc/(iblock-1)]); end;
    
    parasim = zeros(nsim,npara);
    likesim = zeros(nsim,1);
    postsim = zeros(nsim,1);
    rej     = zeros(nsim*ntimes,1);
    TTTsim = zeros(nsim,nstate^2);
    RRRsim = zeros(nsim,nstate*nshocks);
    CCCsim = zeros(nsim,nstate);
    zsim = zeros(nsim,nstate);
    
    for j = 1:nsim*ntimes
        
        Tim = Tim+1;
        
        % Draw para_new from the proposal distribution (a multivariate normal
        % distribution centered on the previous draw, with standard
        % deviation cc*sigscale). 
        
        para_new = para_old + cc*(sigscale*randn(npara,1));
        para_new = para_new.*(1-para_mask)+para_fix.*para_mask;
        [TTT_new,RRR_new,CCC_new,valid_new] = dsgesolv(mspec,para_new,nant);
        retcode = valid_new;
        
        % Solve the model, check that parameters are within bounds, and
        % evalue the posterior. 
        
            [post_new,like_new,zend_new,ZZ_new,DD_new,QQ_new] = feval('objfcnmhdsge_2part',para_new,bounds,YY,YY0,nobs,...
                nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,TTT_new,RRR_new,CCC_new,valid_new,para_mask,coint,cointadd,cointall,YYcoint0,nant,antlags,nshocks-nant);

        
        % Calculate the multivariate log likelihood of jump from para_old to para_new
        propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - 0.5*sigpropdim*log(cc^2) ...
                   -0.5*(para_new - para_old)'*sigpropinv*(para_new - para_old)/cc^2;
        
             
        % Choose to accept or reject the new parameter by calculating the
        % ratio (r) of the new posterior value relative to the old one
        % We compare min(1,r) to a number drawn randomly from a 
        % uniform (0,1) distribution. This allows us to always accept 
        % the new draw if its posterior value is greater than the previous draw's,
        % but it gives some probability to accepting a draw with a smaller posterior value,
        % so that we may explore tails and other local modes.
        
        r = min([1 ; exp( post_new - post_old)]);
        
        if rand(1,1) < r;
            % Accept proposed jump
            para_old = para_new;
            post_old = post_new;
            like_old = like_new;
            
            TTT_old = TTT_new;
            RRR_old = RRR_new;
            CCC_old = CCC_new;
            valid_old = valid_new;
            
            zend_old = zend_new;
            ZZ_old = ZZ_new;
            DD_old = DD_new;
            QQ_old = QQ_new;
            
        else
            % Reject proposed jump
            rej(j) = 1;
            reje = reje+1;
            
        end
        
        % Save every (ntime)th draw
        if (j/ntimes) == round(j/ntimes)
            likesim(j/ntimes,1) = like_old;
            postsim(j/ntimes,1) = post_old;
            parasim(j/ntimes,:) = para_old';
            TTTsim(j/ntimes,:) = TTT_old(:)';
            RRRsim(j/ntimes,:) = RRR_old(:)';
            CCCsim(j/ntimes,:) = CCC_old(:)';
            zsim(j/ntimes,:) = zend_old(:)';
        end
        
        if j == nsim*ntimes
            fprintf(1,'\nRejection perct %2.4f ',[sum(rej)/j]);
        end      
        
    end
     
    fwrite(fid1,parasim','single');
    fwrite(fid2,postsim','single');
    fwrite(fid3,TTTsim','single');
    fwrite(fid4,RRRsim','single');
    fwrite(fid5,zsim','single');
    
end

toc

fprintf(1,'rejection rate = %2.4f',reje/Tim);

%% Calculate Parameter Covariance Matrix

num1 = nblocks*nsim*npara*4;            % number of bites per file
numb1 = nburn*npara*4;                  % number of bites to discard

% read in saved parameter draws
infile_params = [spath,'/params'];
fid = fopen(infile_params,'r');
status = fseek(fid,numb1,'bof');


% Reshape parameter draws into a matrix so can calculate covariance
theta = [];
while ( ftell(fid1) < num1 )
    theta_add = fread(fid,[npara,nsim],'single')';
    theta = [theta; theta_add];
end

% Calculate covariance and save
cov_theta = cov(theta);
fwrite(fid6,cov_theta(:),'single');
fclose('all');
