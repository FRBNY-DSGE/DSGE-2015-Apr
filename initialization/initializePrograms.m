
%% Set Random Number Generators
rand('state',12345)
randn('state',12345)

%% Initialize Model Specifications
spec

%% Load data
[YYall,XXall,ti,nobsall,dlpopall,dlMA_pop,MA_pop,population] = loaddata(nvar,nlags,nant,antlags,psize,zerobound,peachflag,mspec);

% Augment number of observerved variables by number of anticipated policy shocks.
nvar = nvar+nant;

% Define in-sample data
I = find(ti == dates);
YY = YYall(I-stime+1:I,:);
XX = XXall(I-stime+1:I,:);
tiI = ti(I-stime+1:I);
dlpop = dlpopall(I-stime+1:I);
pop_smth = population(I-stime+2+nlags:I+1+nlags);
nobs = size(YY,1);

% Define pre-sample data
YY_p = YYall(1:I-stime,:);
XX_p = XXall(1:I-stime,:);
dlpop_p = dlpopall(1:I-stime);
nobs_p = size(YY_p,1);

% Initialize date-relevant variables
Idate = find(tiI == dates);
tiall    = [tiI(1:end-1);(tiI(end):0.25:(tiI(end)+.25*(qahead+stime-Startdate)))'];
datesall = strcat(num2str(floor(tiall)),'-',num2str(round(1+4*(tiall-floor(tiall)))));

% Transformations associated with level variables
if mspec==557
    YY(1,[1:6 8:end]) = NaN;
    YY(2:end,7)       = NaN;
end

%% Load information for prior distribution
prior = priors990();

pmean  = prior(1:npara,1);
pstdd  = prior(1:npara,2);
pshape = prior(1:npara,3);

pshape = pshape.*(1-para_mask);
pmean = pmean.*(1-para_mask)+para_fix.*para_mask;
pstdd = pstdd.*(1-para_mask);

nonpolipar = [1:1:npara];
nonpolipar(polipar) = [];

%% Load transformation scheme for parameters
trspec = transp990();
trspec(:,1) = trspec(:,1).*(1-para_mask);

%% Set number of states and shocks
[TTT,RRR,CCC,valid] = dsgesolv(mspec, para);
[nstate,nshocks] = size(RRR);

%% Calculate lagged growth rates
YY0 = zeros(nlags,size(YY,2)); 
nobs0 = size(YY0,1);
for lagind = 1:nlags
    YY0(nlags-lagind+1,:) = XX(1,1+cointall+(lagind-1)*nvar+1:1+cointall+lagind*nvar); 
end
YYcoint0=[];
