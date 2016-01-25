% filename: spec.m -- updated with 2010:Q2 data specifications
% description:

%% Main model specifications and important flags
iter = 1;

if ~exist('CH','var'), CH = 1; end
if ~exist('reoptimize','var'), reoptimize = 1; end

%% Number of observations used in VAR
precrisis = 0;


%% Set Number of Anticipated Shocks

% zerobound is a flag for using anticipated shocks to fix interest rate
% expectations near zero. nant is the number of periods that the interest rate
% should be fixed near zero.

if mspec==557
   nant_implied=4;
elseif mspec==990
   nant_implied = 6;
end

if zerobound
    % antlags contains a value for the number of periods back we should start
    % incorporating zero bound expectations.  ZLB expectations should begin
    % 2008Q4
    if ~exist('antlags','var'), antlags = 24; end
    
    % nant contains the number of periods forward that expectations are fixed
    % to zero.
    if ~exist('nant','var'), nant = nant_implied; end
else
    nant = 0;
    antlags = 0;
end

%% Set an interest rate lower bound of 0.25
%The monetary policy shock adjusts to preserve the bound.
bdd_int_rate = 0;

%% Forecast: general options

if ~exist('qahead','var'),qahead = 60; end

%% Forecast: implied observables

history = 0; % plot history (1:stime+qhead)

% alternative specifications: these specs pertain to the model that you
% load actual data for comparison

%% Forecast: conditional data
peachdesclist{1} = 'Unconditional';
peachdesclist{2} = 'Central Scenario';
peachdesclist{3} = 'Cond. FFR & Spread';
peachdesc = char(peachdesclist(peachflag+1));

peachstr = {'';'_peach';'_semicond'};
peachfile=('conddata');

if peachflag
    load(peachfile);
    psize = size(data,1);
    clear data;
else
    psize = 0;
end

%% Parallel Processing

% For alternative parallel programs using vcSubmitQueuedJobs
if ~exist('nMaxWorkers','var'),nMaxWorkers = 20; end
if ~exist('distr','var'); distr = 0; end;
if ~exist('parflag','var'), parflag = 1; end

%% Disturbance Smoothing

%% Turn Shocks off: Forecast distribution then reflects only param uncertainty
if ~exist('sflag','var'), sflag = 0; end % sflag = 1: Turn all shocks off in the forecast

%% Plotting: general specifications (forplot.m and plot_all.m)
if ~exist('useSavedMB','var'), useSavedMB = 0; end
% if useSavedMB = 0, forplot recalculates and resaves means and bands
% if useSavedMB = 1, forplot loads in means and bands saved from a previous run
% means and bands are saved to a file called fcastMeans*.mat.

% when creating this .mat file for the first time, make sure plotList is
% complete so that means and bands are calculated for all figures

if ~exist('fancharts','var'), fancharts = 1; end

% try to get rid of this section
if ~exist('plotSeparate','var'), plotSeparate = 0; end % 1 = each observable is graphed in a separate PDF rather than as part of a 2x2 subplot

plotSeparate_fcast = 1; % 1 = each observable is graphed in a separate PDF rather than as part of a 2x2 subplot
plotSeparate_counter = 0; % 1 = each observable is graphed in a separate PDF rather than as part of a 2x2 subplot
plotSeparate_shockdec = 0;
plotSeparate_irf = 0; % 1 = each observable is graphed in a separate PDF rather than as part of a 2x2 subplot
plotSeparate_shocks = 1; % 1 = each observable is graphed in a separate PDF rather than as part of a 2x2 subplot
plotSeparate_vdec = 1;

%% Plotting: shockdec settings (forplot.m and plot_all.m)
shockdec_history = 0;

%% Plotting: counterfactual settings (forplot.m and plot_all.m)
if ~exist('Counterforecast','var'), Counterforecast = 2; end
if ~exist('Shockremove','var'), Shockremove = 1; end


if ~exist('ShockremoveList','var'),
    switch mspec
        case {557}, ShockremoveList = [0:10];
        case {990}, ShockremoveList = [0:16];
        otherwise, error('Set ShockremoveList variable');
    end
end

if ~exist('stime','var')
    if precrisis==0,
        if mspec== 557
            stime = 122;
        elseif mspec==990
            stime = 220;
        end
    end
end

Enddate_forecastfile = stime+qahead;

%% Plotting: presentation options (forplot.m)
% Note: specifications for make_presentation.m are now in pres_spec.m
if ~exist('plotList','var')
    plotList = {'Forecast';
        'Shock Decomposition';
        'Counterfactual by Variable';
        'Counterfactual by Shock';
        'Shock';
        'ShockAnt'};
    %'Q4Q4 Table';
end

%% Plotting:
percent = 0.68;  % Percent for the bands

if ~exist('fourq_flag','var'), fourq_flag = 0; end


%% mspec_add
% try to incorporate graphing variables into forplot, group with other graphing variables
[nvar,varnames,graph_title,cum_for,popadj,varnames_YL,varnames_irfs,varnames_YL_4Q,varnames_YL_irfs,...
    names_shocks,names_shocks_title,nl_shocks_title,shocksnames,cum_irf,vardec_varnames,shockcats,list,shockdec_color] = mspec_add(mspec,dataset,zerobound,nant,fourq_flag);


cum_for_4q = cum_for;
cum_for_4q(4) = 4;


q_adj=100;


% To plot 99 percent bands
if ~exist('onepctflag','var'), onepctflag = 0; end
if onepctflag
    onepctstr = '1pct';
else
    onepctstr = '';
end

%% Set end date for data (loaddata.m, forplot.m, plot_all.m)
if ~exist('dates','var')
    if precrisis==0, dates = 2014.75; end % for the post-crisis estimate
    if precrisis==1, dates = 2007.25; end % for the pre-crisis estimate
end

if ~exist('mnobss','var'), mnobss = (dates - 2007)*4; end  % use for graphs starting in 2007Q1
pnobss=psize;
plot_forward = 13;

future = plot_forward+4-round(1+4*(dates-floor(dates)));

if exist('shockdec_history','var') && shockdec_history==1
    Startdate = 1;
else
    Startdate= stime - mnobss;
end

if Counterforecast==0
    Enddate = stime + pnobss;
elseif Counterforecast==1
    Enddate = stime + qahead;
elseif Counterforecast==2 % Enddate gets set in forplot_dick % Moving it back to spec to avoid duplication of vars in make_pres
    Enddate =stime + future;
end

counter_ahead = (Enddate_forecastfile - Startdate);

%% Set start and end dates for X-axis
% sirf different for each type of figure
% sirf needs to be created here so that means and bands have correct pre-allocated size
% remove sirf from figspecs
sirf = (Startdate:stime+future); % Startdate and future are set in spec
sirf_shockdec = (Startdate:Enddate);
sirf_counter = (Startdate:Enddate);
sirf_shock_1 = (Startdate:stime); % further modified in plotNL

%% IRF Settings (irfsim.m, irfplot.m)
% redo_irfsims similar to useSavedMB for parforecast.
% 1 = re-calculate IRFs using draws from gibb
% 0 = re-caculate means and bands using existing IRFs
if ~exist('redo_irfsims','var'), redo_irfsims = 1; end
if ~exist('nirf','var'), nirf = 40; end
nplotstates = 18;
if ~exist('irfStates','var'), irfStates = 0; end
fix = -0.5; % Calculate a <fix> bp policy shock. If fix is negative, then this calculates an expansionary policy shock.

%% Moments options (mom.m)
paragroupseq = [''];
if ~exist('PLOTDRAWS','var'), PLOTDRAWS = 0; end
%paragroupseq = ['' 'l' 'g']; % to see moments for only those parameters
% for which the inflation forecast 20 quarters out is
% greater than 4 (g) or less than 3 (l)

%% Configure the Metropolis Algorithm (number of simulations) (gibb.m, parforecast.m, forplot.m)
if ~exist('nblocks','var'), nblocks = 11; end
if ~exist('nsim','var'), nsim = 10000; end
if ~exist('nburn','var'), nburn = nsim; end % initial simulations to 'burn'
if ~exist('jstep','var'), jstep = 5; end  % decreasing jstep to make forecasts and bands smoother
ntimes = 5;

%% mspec_parameters; mspec_parameters_'mspec'(subspec,dataset);
eval(['[para,para_names,para_mask,para_fix,npara,polipar,polivalue,bounds] = mspec_parameters_',num2str(mspec),'(',num2str(subspec),',',num2str(dataset),');']);

%% Name string variables: Can be used in saving
lmodel = ['m',num2str(mspec),num2str(subspec)];
lprior = [num2str(mprior),num2str(pf_mod)];
ds = num2str(10*dataset);
if sflag==0;
    ssf='';
else
    ssf=num2str(sflag);
end

if isempty(stime),
    st = '';
else
    st = num2str(stime);
end
if mspec == 1, st = num2str([]); end;

if zerobound
    antstr = ['ZB_' num2str(nant) '_L' num2str(antlags)];
    antstrold = ['ZB_' num2str(nant+1) '_L' num2str(antlags-1)];
else
    antstr = '';
    antstrold = '';
end

if parflag, parstr = 'par'; else parstr = ''; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFREQUENTLY CHANGED VARIABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specific to Euro, we don't understand

coint = 0;
cointadd = 0;
cointall = coint + cointadd;

%% choose policy option
Policy_Sel = [];

%% use prior (= 1) or posterior (= 0)
if ~exist('PRIO','var')
    PRIO = 0;
end
%% Flag - exclude monetary policy shocks from loss computation
NoMon = 1;

%% Flag - allow for mispecification in the estimated policy rule (backward)
MPol = [];
if MPol
    MPol_ = 'Mpol';
else
    MPol_ = '';
end

%% change deep parameters, otherwise set Ideep = []
Ideep = [];%[9 .6];

%% number of lags
if mspec==557
    nlags = 0;
elseif mspec == 990
    nlags = 2;
else
    error('nlags needs to be defined for specified mspec.');
end

%% End of Presample
%T0 = nlags; %%0
T0 = 0; % Set this to nlags unless you want no presample - then do 0.
Tend = 0; % End of Estimation Sample

%% note: max number of inderminancies across all draws/policy parameter combinations.
%% if actual > maxind for some draws/policy parameter combination , programs won't work
%% if storage space is not concern, just set maxind high (say 5 or so)
IND = 1;
maxind = 3;


