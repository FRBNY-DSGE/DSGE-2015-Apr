%% MODEL-SPECIFIC SPECIFICATIONS (Defaults)
  mspec   = 990;  % Model to run
  subspec = 2;    % Subspecifications to change parameter values
  mprior  = 9;    % Prior Specifications
  pf_mod  = '90'; % Prior Specifications
  dataset = 715;  % Dataset number; can be used to reference alternative
                  % datasets
  nant    = 6;    % Number of anticipated policy shocks

%% SAVE SETTINGS
  overwrite = 1; % Overwrite old output?

%% ESTIMATION SETTINGS
  reoptimize = 0;     % Whether to reoptimize or use the mode as is
  CH         = 0;     % Compute Hessian again
  nsim       = 10000; % Number of draws per posterior simulation block
  nblocks    = 22;    % Number of blocks
  nburn      = nsim*2;  % How many draws to discard as burn-in
  jstep      = 10;     % From the blocks, we take every jstep-th element

%% FORECAST SETTINGS

  zerobound = 1; % Incorporate anticipated policy shocks
  antpolflag=1;

  % CONDITIONAL DATA
  %  Indicates whether we are doing a forecast conditional on an
  %  additional quarter of data for certain variables.
  peachflag = 1;

  % PARALLEL SETTINGS
  distr = 1;      % Run forecasting in parallel; requires parallel toobox
  nMaxWorkers=40; % If run in parallel, how many workers?

