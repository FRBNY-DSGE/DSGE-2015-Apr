%                                   OVERVIEW
% forecast_parallel.m: Computes various forecasts using vcSubmitQueuedJobs.m to
%                      distribute draws among workers. Each worker
%                      constructs forecasts using forecastFcn.m.
%
%                      This program uses draws from gibb.m to create forecasts,
%                      counterfactuals,shock decompositions, and smoothed
%                      shock
%                      estimates, both unconditional and conditional on
%                      peach Peach's forecast.
%
%                           IMPORTANT SUBFUNCTIONS
% setInfiles.m, setOutfiles.m: specify infiles and outfiles
% forecastFcn.m: subfunction distributed to each worker, to compute
%                forecasts.

%% Inititalization
tic;
close all;
initializePrograms;

%% Set warning off for newer versions of MATLAB
try
  warning('off', 'MATLAB:warn_r14_stucture_assignment')
end

%% Specifications

% Set dimensions of pre-allocated matrices
nstate_0  = nstate;
nshocks_0 = nshocks;

% Augment the list of shocks to add the anticipated policy shocks
if zerobound == 1
    nshocks_all = nshocks+nant;
    ShockremoveList = [ShockremoveList,(nshocks+1:nshocks+nant)];
else
    nshocks_all = nshocks;
end

%% Read infiles (outputted from gibb_est_ant.m).

st_adj = adj2part(mspec);

% List data to be read in
indataType = {'params','TTT','RRR','zend'};

% Create structures infile (filenames), numb (# bytes to burn) , and num (# bytes per file)
[infile,numb,num] = setInfiles(indataType,...
    spath,nburn,nstate+nant_implied+st_adj,nshocks+nant_implied,npara,nblocks,nsim);
listInfiles = fieldnames(infile);

for iInfile = 1:length(listInfiles)
    infileName = listInfiles{iInfile};
    fid.(infileName) = fopen(infile.(infileName),'r');
    status.(infileName) = fseek(fid.(infileName),numb.(infileName),'bof');

end
lambdas=NaN;

st_adj = adj2part(mspec);

%% Open outfiles

% Use spath_overwrite if you want to save output files to a different
% folder than spath.
if exist('spath_overwrite','var'), spath = spath_overwrite; end

% List data to be written
outdataType = {'hist','forecast','forecast_s'};
if peachflag,
    outdataType = union(outdataType,{'condhist','condforecast','condforecast_s'});
end
outdataType = union(outdataType,{'states','condforecast','peachshocks','condshocks','datashocks',...
    'condforecast_ns','peachshocks_ns','condshocks_ns','datashocks_ns',...
    'counter','shockdec','ytrend','dettrend','dettrend_peach'});

% Create outfile structure (output filenames)
[outfile] = setOutfiles(outdataType,zerobound, spath,peachflag,ShockremoveList);

% Open outfiles to write
if overwrite
    listOutfiles = fieldnames(outfile);
    for iOutfile = 1:length(listOutfiles)
        outfileName = listOutfiles{iOutfile};
        fid.(outfileName) = fopen(outfile.(outfileName),'w');
    end
end

%% Load in peachdata and specify qahead_extra
% peachdata (nant x (nvar+nant)): a matrix of conditional data (usually peach Peach's
%                                 nowcast). If we are doing a conditional forecast,
%                                 we append peachdata to the end of our observed data.
% qahead_extra: If we are doing a conditional forecast, we forecast
%               qahead_extra (59Q) quarters ahead (psize (1Q) quarters less
%               than the unconditional forecast) instead of qahead (60Q).

% load conditional data
if peachflag > 0,
    load(peachfile);
    peachdata = data;
    
    % Population adjust and get into logs
    peachdata(1) = peachdata(1)-400*dlMA_pop(1);
    
    %for models which read in quarterly data, adjust peachdata (Except
    %hours)
    peachdata([1 3:end]) = peachdata([1 3:end])/4;
    
    if size(peachdata,2) < nvar
        peachdata = [peachdata,NaN(size(peachdata,1),nvar-size(peachdata,2))];
    end
    
    if zerobound

        [ExpFFR,peachdata_FFR] = ExpFFR_OIS(nant,antlags,psize,zerobound,peachflag);

        ExpFFR = ExpFFR/4;
        peachdata_FFR = peachdata_FFR/4;

        if ~is2part(mspec)
            peachdata = [peachdata,NaN(size(peachdata,1),nant)];
        end

        for t = 1:size(peachdata,1)
            peachdata(t,end-nant+1:end-t) = peachdata_FFR(t,1:nant-t);
        end
    else
        ExpFFR = [];
    end
    
    qahead_extra = qahead - size(peachdata,1);
else
    if zerobound
        [ExpFFR,peachdata_FFR] = ExpFFR_OIS(nant,antlags,psize,zerobound,peachflag);
        if class_mspec(mspec) == 904 %for models which read in quarterly data, adjust FFRs
            ExpFFR = ExpFFR/4;
            peachdata_FFR = peachdata_FFR/4;
        end
    else
        ExpFFR = [];
    end
    peachdata = [];
    qahead_extra = qahead;
end

semi_peachdata = [];


%% Load in parameter draws, state transition equation matrices, and zend (output from gibb.m)
% priotheta (npara x nsim): parameter draws
% TTTsim (nstate^2 x nsim), RRRsim (nstate*nshocks x nsim), CCCsim (nstate x nsim): matrices of the state transition equation S_t = TTT*S_(t-1) + RRR*eps_t + CCC
% zendsim (nstate x nsim): end-of-sample smoothed state: S_(T/T)

% priotheta, TTTsim, RRRsim, CCCsim, and zendsim are compressed matrices.
% They are read in by block (with nsim simulations per block) and reshaped.

priotheta = [];
iblock = 0;
% Read blocks of size nsim
while ( ftell(fid.params) < num.params ) % only if PRIO==0
    
    nstate = nstate_0;
    nshocks = nshocks_0;
    
    iblock = iblock+1;
    fprintf(1,' block: %2.0f; \n',iblock);
    
    % only if PRIO==0
    priotheta = fread(fid.params,[npara,nsim],'single')';
    TTTsim = fread(fid.TTT,[(nstate+nant_implied+st_adj)^2,nsim],'single')';
    RRRsim = fread(fid.RRR,[(nstate+nant_implied+st_adj)*(nshocks+nant_implied),nsim],'single')';
    zendsim = fread(fid.zend,[(nstate+nant_implied+st_adj),nsim],'single')';
    
    if ~zerobound
        priotheta1 = priotheta(1:end-20,:);
        
        nstate_sim = nstate+nant_implied+1;
        nshocks_sim = nshocks+nant_implied;
        iSplit = getState(mspec,0,'n_end+n_exo');
        
        selectorTTT = ones(nstate_sim, nstate_sim);
        selectorTTT(:,iSplit+1:iSplit+(nant_implied+1)) = 0;
        selectorTTT(iSplit+1:iSplit+(nant_implied+1),:) = 0;
        selectorTTT = selectorTTT(:);
        TTTsim = TTTsim(:,logical(selectorTTT));
        
        selectorRRR = ones(nstate_sim, nshocks_sim);
        selectorRRR(iSplit+1:iSplit+(nant_implied+1),:) = 0;
        selectorRRR(:,end-nant_implied+1:end) = 0;
        selectorRRR = selectorRRR(:);
        RRRsim = RRRsim(:,logical(selectorRRR));
        
        selectorzend = ones(nstate_sim,1);
        selectorzend(iSplit+1:iSplit+(nant_implied+1)) = 0;
        zendsim = zendsim(:,logical(selectorzend));
    end
    
    priotheta = priotheta(jstep:jstep:nsim,:);

    TTTsim = TTTsim(jstep:jstep:nsim,:);
    RRRsim = RRRsim(jstep:jstep:nsim,:);
    zendsim = zendsim(jstep:jstep:nsim,:);

    %% Distribute draws to nMaxWorkers, as input for forecastFcn.m
    % idx: how many draws to distribute to each worker
    % JobOptions: holds input (to function irfsimFcn.m) for each of the nMaxWorkers
    % nout: number of output arguments
    % pathdep.m: a program that establishes variable PathDependencies
    % JobOut: holds output from forecastFcn.m
    valid_para=1;
    
    if distr
        nout = 6;
        
        JobOptions = cell(1,nMaxWorkers);
        for k = 1:nMaxWorkers
            idx = ceil(nsim/(jstep*nMaxWorkers))*(k-1)+(1:ceil(nsim/(jstep*nMaxWorkers)));
            idx = idx(idx<=nsim/jstep);
            JobOptions{k} = {k,priotheta(idx,:),TTTsim(idx,:),RRRsim(idx,:),zendsim(idx,:),...
                mspec,zerobound,bdd_int_rate,peachflag,0,0,...
                nvar,qahead,qahead_extra,nshocks_all,YY,nobs,peachdata,psize,...
                Enddate_forecastfile,Startdate,stime,ShockremoveList,...
                nant,antlags,varnames,names_shocks,...
                nstate_0,nshocks_0,sflag,YY0,...
                valid_para,nlags,npara,coint,cointadd,ExpFFR,nplotstates};
            
            pathdep;
        end
        
        
        JobOut = vcSubmitQueuedJobs(nMaxWorkers,@forecastFcn_est_ant,...
          nout,JobOptions,PathDependencies,nMaxWorkers);
        
        % Concatenate output across workers
        temp_fcast = [JobOut{:,1}];
        fcast.uncond = cat(1,temp_fcast(:).uncond);
        if zerobound, fcast.uncond_bdd = cat(1,temp_fcast(:).uncond_bdd); end
        if peachflag
            fcast.cond = cat(1,temp_fcast(:).cond);
            if zerobound, fcast.cond_bdd = cat(1,temp_fcast(:).cond_bdd); end
        end
            
        temp_shocks = [JobOut{:,2}];
        shocks.data = cat(1,temp_shocks(:).data);
        shocks.data_ns = cat(1,temp_shocks(:).data_ns);
        if peachflag
            shocks.peach = cat(1,temp_shocks(:).peach);
            shocks.datawpeach = cat(1,temp_shocks(:).datawpeach);
            shocks.peach_ns = cat(1,temp_shocks(:).peach_ns);
            shocks.datawpeach_ns = cat(1,temp_shocks(:).datawpeach_ns);
        end
        ytrend = cat(1,JobOut{:,3});
        dettrend = cat(1,JobOut{:,4});
        counter_all = cat(1,JobOut{:,5});
        shockdec_all = cat(1,JobOut{:,6});
        
        clear JobOut JobOptions;
        
    else
        
        % Compute forecasts for each draw sequentially (if debugging, for example)
        k=1;
        [fcast,shocks,ytrend,dettrend,counter_all,shockdec_all] = ...
          forecastFcn_est_ant(k,priotheta,TTTsim,RRRsim,zendsim,...
            mspec,zerobound,bdd_int_rate,peachflag,0,0,...
            nvar,qahead,qahead_extra,nshocks_all,YY,nobs,peachdata,psize,...
            Enddate_forecastfile,Startdate,stime,ShockremoveList,...
            nant,antlags,varnames,names_shocks,...
            nstate_0,nshocks_0,sflag,YY0,...
            valid_para,nlags,npara,coint,cointadd,ExpFFR,nplotstates);
        
    end
    
    %% Save forecasts, and other relevant output
    if overwrite
        fwrite(fid.forecast,fcast.uncond','single');
        if zerobound, fwrite(fid.forecastbdd,fcast.uncond_bdd','single'); end
        fwrite(fid.datashocks,shocks.data','single');
        fwrite(fid.datashocks_ns,shocks.data_ns','single');
        for peachcounter = 0:peachflag
          shockcount = 0;
          for Shockremove = ShockremoveList
            shockcount = shockcount + 1;
            fwrite(fid.(['counter',num2str(peachcounter),num2str(Shockremove)]),counter_all(:,:,peachcounter+1,shockcount)','single');
          end
        end
        for peachcounter = 0:peachflag
          shockcount = 0;
          for Shockremove = ShockremoveList
            if Shockremove > 0
              shockcount = shockcount + 1;
              fwrite(fid.(['shockdec',num2str(peachcounter),num2str(Shockremove)]),shockdec_all(:,:,peachcounter+1,shockcount)','single');
            end
          end
        end
        
        if peachflag
          fwrite(fid.peachshocks,shocks.peach','single');
          fwrite(fid.condshocks,shocks.datawpeach','single');
          fwrite(fid.peachshocks_ns,shocks.peach_ns','single');
          fwrite(fid.condshocks_ns,shocks.datawpeach_ns','single');
          fwrite(fid.dettrend_peach,dettrend(:,:,2)','single');
          fwrite(fid.condforecast,fcast.cond','single');
          if zerobound, fwrite(fid.condforecastbdd,fcast.cond_bdd','single'); end
        end
        fwrite(fid.ytrend,ytrend','single');
        fwrite(fid.dettrend,dettrend(:,:,1)','single');
    end
end

fclose('all');

fprintf('\n Elapsed time is %4.2f minutes',toc/60);
