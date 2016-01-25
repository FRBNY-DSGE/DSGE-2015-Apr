% OVERVIEW
% keepVars.m: Keep important variables in the workspace. This allows us to
%             loop over programs. 
%
% IMPORTANT VARIABLES
% keepList: a list of variables to keep

%% Important
keepList = {'overwrite', 'reoptimize'};

%% Model 
keepList = union(keepList,{'mspec','subspec','pf_mod','mprior','peachsm','dataset','dates','stime', ...
                            'altdata'});

%% Zerobound 
keepList = union(keepList,{'zerobound','nant','nant_implied','antlags'});

%% Parallel
keepList = union(keepList,{'parallelflag','parflag','distr','nMaxWorkers'});

%% Conditional 
keepList = union(keepList,{'peachfile','peachfileold'});

%% Forecasting 
keepList = union(keepList,{'simple_forecast','bdd_int_rate','peachflag','dsflag','jstep','save_states','qahead','nplotstates'});

%% Forplot
keepList = union(keepList,{'plot*','fancharts','fourq_flag','Q4Q4table','q_adj', 'vars_to_plot'});

%% Plotting 
keepList = union(keepList,{'newsletter','system','pres','useSavedMB','plotList','medianFlag','issue_num'});

%% IRFs
keepList = union(keepList,{'irfStates','redo_irfsims','nirf'});

%% Vardecs
keepList = union(keepList,{'vdec_states'});

%% PLT
keepList = union(keepList,{'MODEstr','experiment_flag','useMode'});

%% Debugging
keepList = union(keepList,{'nsim','nblocks','nburn'});

%% Paths
keepList = union(keepList,{'spath','fpath','spath_overwrite','fpath_overwrite','gpath_overwrite'});

%% Output-producing scripts
keepList = union(keepList,{'baseDir', 'matlabDir', 'spec_file','systemDir', 'my_keepList', 'forecast_again'});

%% For script-wrapping
keepList = union(keepList,{'varargin', 'script'});
    
%% Productivity
keepList = union(keepList,{'prodflag'});

%% Aux variables
keepList = union(keepList,{'auxFlag','bolUseFinals','mnobss','shockdec_history','useRtParallelResults','saveSemiCondStShFlag',...
    'paramFile','paramSpecFile', 'replaceBlock','edate','sdate','h1','h2','h3','hairType','haircompare','swapKapSSS','rwMC','varName',...
    'sh_ind','st_ind','sh_cmp','st_cmp'});

%% Misc
% iter save_states iInput stime dates mnobss shockdec_history
% experiment_flag useMode MODEstr input
% irfStates vdec_states sd_states
% modeonly ssonly ShockremoveList

%% T-distributed shocks
keepList = union(keepList,{'tflag','df_bar_num','dfflag','pargibbflag','onepctflag','tstr'});

%% Realtime 
keepList = union(keepList,{'judge','qvint','ovint','sflag','r_exp','nantmax','antpolflag','mant','first','compound','paramDate'});

%% PLT specs
if exist('PLT','var') && PLT==1
    keepList = union(keepList,{'experiment_flag','useMode','nsim','nburn','nblocks',...
                               'jstep','init*','pigap_0','MODEstr','EXPE','exp','expers',...
                               'plottype','fancharts','Means_*','alt_baseline'});
end

%% Forecast Loop through multiple models (see forecastAndPlotLoop.m)
keepList = union(keepList,{'cellSpecList'});

%% Alternate Gensys
keepList = union(keepList,{'gensys2', 'model2_num', 'add', 'mspec1', 'mspec2', 'infile0', 'param_555', 'param_556'});

%% IRF one draw and output gap
keepList = union(keepList,{'irfvars','irfstates','ygapflag','params'});

%% Forward guidance dates, colors, plot data
keepList = union(keepList,{'exercise_flag','policy_flag','t1','t2','t2path','t2i','FGplot','FGbaseline'});
keepList = union(keepList,{'iPlot','types','type','colors','ind','inds'});

% alternative policy rules
keepList = union(keepList,{'alt_policy', 'alt_policy_flag', 'alt_policy_list', 'alt_policy_series', 'rule', 'ii', 'altrules'});
keepList = union(keepList, {'alt_policy_spath_override'});
keepList = union(keepList, {'to_save', 'start_loc', 'quarters_fwd', 'folder_loc'});

% specific plot-path
keepList = union(keepList, {'gpath_override'});
keepList = union(keepList, {'gpath', 'path_to_plotdirs'});
keepList = union(keepList, {'scenario_fcast'});

% I want to use variables that don't disappear!! - MJS 08-01-14
keepList = union(keepList, {'myCustomVar1', 'myCustomVar2', 'myCustomVar3'});

eval(['keep ', sprintf(repmat('%s ',1,length(keepList)),keepList{:})])
