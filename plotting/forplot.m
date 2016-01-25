% OVERVIEW (written Jan 23, 2014 by MDC, day this is put into use)
% forplot.m: Computes the means and bands for the states; prepares
%                   for plotting.
% IMPORTANT VARIABLES
% Means: a structure variable with fields for each type of plot. Each field
%        holds the mean <forecast/counterfactual/etc> across draws.
% Bands: a structure variable with fields for each type of plot. Each field
%        holds fanchart bands across draws.
% useSavedMB: a flag that allows us to use means and bands calculated and
% saved from a
%             previous run of forplot. This option allows us to avoid
%             re-loading in forecasts computed across draws, which is the most
%             time-consuming aspect of forplot.
% plotList: a cell array containing the types of plots you want to
%           produce. If not set, this will be set to a default by spec.m 


clear infile*
close all;
keepVars;
initializePrograms;


% Set up 'number of' dimension to either nvar or nplotstates
% Set up the suffix that imports either the states or the observable forecast

    ndim = nvar;
    
    cum_for_plot = cum_for;
    popadj_plot = popadj;


notzero = ShockremoveList > 0;
ShockremoveList = ShockremoveList(notzero);

if zerobound
    ShockremoveList = [ShockremoveList,(nshocks+1:nshocks+nant)];
    nshocks = nshocks+nant;
end


record = struct([]);

%% Re-calculate means and bands
if ~useSavedMB
    
    %% Initialize means and bands
    Means.forecast = zeros(ndim,qahead,peachflag+1);
    Bands.forecast = zeros(2,ndim,qahead,peachflag+1);
    
    if onepctflag, Bands_onepct = zeros(2,ndim,qahead,peachflag+1); end
    
    if any(ismember(plotList,{'Counterfactual by Variable', 'Counterfactual by Shock'}))
        Means.counter = zeros(ndim,counter_ahead,peachflag+1,length(ShockremoveList));
        Bands.counter = zeros(2,ndim,counter_ahead,peachflag+1,length(ShockremoveList));
    end
    
    if any(ismember(plotList,{'Shock Decomposition','Shock Decomposition-noDet'}))
        Means.trend = zeros(ndim,counter_ahead);
        Bands.trend = zeros(2,ndim,counter_ahead);
        
        Means.shockdec = zeros(ndim,counter_ahead,peachflag+1,length(ShockremoveList));
        Bands.shockdec = zeros(2,ndim,counter_ahead,peachflag+1,length(ShockremoveList));
        
        Means.dettrend = zeros(ndim,counter_ahead,peachflag+1);
    end
    
    if any(ismember(plotList,{'Shock'}))
        Means.uncondshocks = zeros(nshocks,nobs);
        Bands.uncondshocks = zeros(2,nshocks,nobs);
        if peachflag
            Means.condshocks = zeros(nshocks,nobs+psize);
            Bands.condshocks = zeros(2,nshocks,nobs+psize);
        end
    end
    
    if any(ismember(plotList,{'ShockAnt'}))
        Means.uncondshocks_ns = zeros(nshocks,nobs);
        Bands.uncondshocks_ns = zeros(2,nshocks,nobs);
        if peachflag
            Means.condshocks_ns = zeros(nshocks,nobs+psize);
            Bands.condshocks_ns = zeros(2,nshocks,nobs+psize);
        end
    end
    
    if any(ismember(plotList,{'4Q Forecast'}))
        Means.forecast4q = zeros(size(Means.forecast));
        Bands.forecast4q = zeros(size(Bands.forecast));
    end
    
    if any(ismember(plotList,'Q4Q4 Table'))
        Means.table4q = zeros(size(Means.forecast));
        Bands.table4q = zeros(size(Bands.forecast));
        
        Means.tableq4q4 = zeros(ndim,3,peachflag+1);
        Bands.tableq4q4 = zeros(2,ndim,3,peachflag+1);
    end
    
    %% Open forecast files to read from
    
    if any(ismember(plotList,{'Forecast','4Q Forecast','Counterfactual by Variable','Counterfactual by Shock','Shock Decomposition','Shock Decomposition-noDet','Q4Q4 Table'}))
        if zerobound
            infile6 = [spath,'/forecast'];
            infile16 = [spath,'/forecastbdd'];
            if peachflag
                infile1 = [spath,'/condforecast'];
                infile11 = [spath,'/condforecastbdd'];
            end
        else
            infile6 = [spath,'/forecast'];
            if peachflag
                infile1 = [spath,'/condforecast'];
            end
        end
    end
    
    if any(ismember(plotList,{'Shock'}))
        infile4 = [spath,'/datashocks'];
        if peachflag
            infile2 = [spath,'/peachshocks'];
            infile3 = [spath,'/condshocks'];
        end
    end
    
    if any(ismember(plotList,{'ShockAnt'}))
        infile45 = [spath,'/datashocks_ns'];
        if peachflag
            infile25 = [spath,'/peachshocks_ns'];
            infile35 = [spath,'/condshocks_ns'];
        end
    end


    if any(ismember(plotList,{'Counterfactual by Variable', 'Counterfactual by Shock'}))
      for peachcount = 0:peachflag
        for Shockremove = ShockremoveList
          if peachcount
              eval(['infile1' num2str(peachcount) num2str(Shockremove) ' = [spath,''/counter_peach'',' 'num2str(Shockremove)' '];'])
          else
              eval(['infile1' num2str(peachcount) num2str(Shockremove) ' = [spath,''/counter'',' 'num2str(Shockremove)' '];'])
          end
        end
      end
    end
    

    if any(ismember(plotList,{'Shock Decomposition','Shock Decomposition-noDet'}))
        infile200 = [spath,'/ytrend'];
        for peachcount = 0:peachflag
          for Shockremove = ShockremoveList
            if peachflag
              eval(['infile2' num2str(peachcount) num2str(Shockremove) ' = [spath,''/shockdec_peach'', ' 'num2str(Shockremove)' '];']);
            else
              eval(['infile2' num2str(peachcount) num2str(Shockremove) ' = [spath,''/shockdec'',' 'num2str(Shockremove)' '];']);
            end
          end
        end
        infile300 = [spath,'/dettrend'];
        infile301 = [spath,'/dettrend_peach'];
        infile301_semi = [spath,'/dettrend_semipeach'];
    end
    
    
    %% Specify num_* (number of bytes per file) and numb_*(number of bytes to discard)
    
    numb_for = nburn*(ndim)*qahead*4;
    num_for = (nblocks*nsim*(ndim)*qahead*4-numb_for)/jstep;
    
    numb_counter = nburn*ndim*counter_ahead*4;
    num_counter = (nblocks*nsim*ndim*counter_ahead*4-numb_counter)/jstep;
    
    numb_shockdec = nburn*ndim*(counter_ahead+1)*4;
    num_shockdec = (nblocks*nsim*ndim*(counter_ahead+1)*4-numb_shockdec)/jstep;
    
    numb_trend = nburn*ndim*4;
    num_trend = (nblocks*nsim*ndim*4 - numb_trend)/jstep;
    
    numb_forshocks = nburn*nshocks*psize*4;
    num_forshocks = (nblocks*nsim*nshocks*psize*4-numb_forshocks)/jstep;
    
    numb_data = nburn*nshocks*nobs*4;
    num_data = (nblocks*nsim*nshocks*nobs*4-numb_data)/jstep;
    
    numb_hist = nburn*(ndim)*stime*4;
    num_hist = (nblocks*nsim*(ndim)*stime*4-numb_hist)/jstep;
    
    %% Open infiles to read
    % infile should not be a variable name, if used for a purpose other than this
    % fid<num> should correspond exactly to infile<num>,
    % ie fid1=fopen(infile1,'r'), not fid01 = fopen(infile1,'r')
    listInfiles = who('infile*');
    for nInfiles = 1:size(listInfiles)
        fidNum = strrep(listInfiles{nInfiles},'infile','');
        eval(['fid',fidNum,'=fopen(infile',fidNum,',','''r'');'])
    end
    
    %% Population adjustment
    % We add Macroeconomic Advisers' population growth forecast to
    % transform our forecasts from per-capita to aggregate. For forecast
    % horizons beyond the scope of MA's forecasts, we use MA's furthest
    % horizon forecast.
    if any(popadj == 1)
        if any(ismember(plotList,{'Forecast','4Q Forecast','Q4Q4 Table'}))
            if qahead > size(dlMA_pop,1)
                dlpop_fcst = [repmat(dlMA_pop',((nblocks-nburn/nsim)*nsim/jstep),1), repmat(dlMA_pop(end),((nblocks-nburn/nsim)*nsim)/jstep, qahead - size(dlMA_pop,1))];
            else
                dlpop_fcst = repmat(dlMA_pop(1:qahead)',(nblocks*nsim)/jstep,1);
            end
        else
            dlpop_fcst = [];
        end
        
        if any(ismember(plotList,{'4Q Forecast','Q4Q4 Table','Counterfactual by Variable','Counterfactual by Shock','Shock Decomposition','Shock Decomposition-noDet'}))
            if counter_ahead - mnobss > size(dlMA_pop,1)
                dlpop_ctr = [repmat(dlMA_pop',((nblocks-nburn/nsim)*nsim/jstep),1), repmat(dlMA_pop(end),((nblocks-nburn/nsim)*nsim)/jstep,counter_ahead-mnobss-size(dlMA_pop,1))];
            else
                dlpop_ctr = repmat(dlMA_pop(1:counter_ahead-mnobss)',(nblocks*nsim)/jstep,1);
            end
            dlpop_old = repmat(dlpop(end-mnobss+1:end)',((nblocks-nburn/nsim)*nsim/jstep),1);
        else
            dlpop_ctr = [];
        end
    else
        dlpop_fcst=[];
        dlpop_old=[];
        dlpop_ctr = [];
    end
    
    % Read in implied state histories if we're plotting states
    
    
    %% Forecast
    
    % Read in data
    if any(ismember(plotList,{'Forecast','4Q Forecast','Counterfactual by Variable','Counterfactual by Shock','Shock Decomposition','Shock Decomposition-noDet','Q4Q4 Table'}))
        ypath = [];
        while ( ftell(fid6) < num_for ) % Read blocks of size nsim
            ypathadd = fread(fid6,[(ndim)*qahead,nsim/jstep],'single')';
            ypath = [ypath;ypathadd];
            ftell(fid6);
        end;
        clear ypathadd;
        fclose(fid6);
        
        if zerobound == 1
            ypath_bdd = [];
            while (ftell(fid16) < num_for)
                ypathadd = fread(fid16,[(ndim)*qahead,nsim/jstep],'single')';
                ypath_bdd = [ypath_bdd;ypathadd];
                ftell(fid16);
            end
            clear ypathadd
            fclose(fid16);
        end
        disp('Loaded Unconditional Forecast')
    end
    
    if any(ismember(plotList,{'Forecast','4Q Forecast','Q4Q4 Table'})) && peachflag
        ypath_cond = [];
        while ( ftell(fid1) < num_for ) % Read blocks of size nsim
            ypathadd = fread(fid1,[(ndim)*qahead,nsim/jstep],'single')';
            ypath_cond = [ypath_cond;ypathadd];
            ftell(fid1);
        end;
        clear ypathadd;
        fclose(fid1);
        
        if zerobound == 1
            ypath_cond_bdd = [];
            while (ftell(fid11) < num_for)
                ypathadd = fread(fid11,[(ndim)*qahead,nsim/jstep],'single')';
                ypath_cond_bdd = [ypath_cond_bdd;ypathadd];
                ftell(fid11);
            end
            clear ypathadd
            fclose(fid11);
        end
        disp('Loaded Conditional Forecast')
        
    end
    
    % Calculate means and bands
    
    if any(ismember(plotList,{'Forecast','4Q Forecast','Counterfactual by Variable','Counterfactual by Shock','Shock Decomposition','Shock Decomposition-noDet','Q4Q4 Table'}))
        for peachcount = 0:peachflag
            for Ivar = 1:ndim
                
                switch peachcount
                    case 1
                        yp = getyp(ypath_cond,qahead,Ivar,cum_for_plot,popadj_plot,dlpop_fcst,YY,Idate,q_adj);
                    case 0
                        yp = getyp(ypath,qahead,Ivar,cum_for_plot,popadj_plot,dlpop_fcst,YY,Idate,q_adj);
                end
                
                if exist('medianFlag','var') && medianFlag==1
                    Means.forecast(Ivar,:,peachcount+1) = median(yp,1);
                else
                    Means.forecast(Ivar,:,peachcount+1) = mean(yp,1);
                end
                
                
                if fancharts
                    Bands.forecast([1 10],Ivar,:,peachcount+1) = hpdint(yp,0.90,0);
                    Bands.forecast([2 9],Ivar,:,peachcount+1) = hpdint(yp,0.80,0);
                    Bands.forecast([3 8],Ivar,:,peachcount+1) = hpdint(yp,0.70,0);
                    Bands.forecast([4 7],Ivar,:,peachcount+1) = hpdint(yp,0.60,0);
                    Bands.forecast([5 6],Ivar,:,peachcount+1) = hpdint(yp,0.50,0);
                else
                    Bands.forecast(:,Ivar,:,peachcount+1) = hpdint(yp,percent,0);
                end
                
                if onepctflag, Bands_onepct(:,Ivar,:,peachcount+1) = hpdint(yp,0.99,0); end
                
            end
            
            if zerobound
                for Ivar = 1:ndim
                    switch peachcount
                        case 1
                            yp = getyp(ypath_cond_bdd,qahead,Ivar,cum_for_plot,popadj_plot,dlpop_fcst,YY,Idate,q_adj);
                        case 0
                            yp = getyp(ypath_bdd,qahead,Ivar,cum_for_plot,popadj_plot,dlpop_fcst,YY,Idate,q_adj);
                            
                    end
                    if fancharts
                        Bands.forecast([1 10],Ivar,:,peachcount+1) = hpdint(yp,0.90,0);
                        Bands.forecast([2 9],Ivar,:,peachcount+1) = hpdint(yp,0.80,0);
                        Bands.forecast([3 8],Ivar,:,peachcount+1) = hpdint(yp,0.70,0);
                        Bands.forecast([4 7],Ivar,:,peachcount+1) = hpdint(yp,0.60,0);
                        Bands.forecast([5 6],Ivar,:,peachcount+1) = hpdint(yp,0.50,0);
                    else
                        Bands.forecast(:,Ivar,:,peachcount+1) = hpdint(yp,percent,0);
                    end
                end
                
                ind_r = find(strcmp(varnames,'Interest Rate'));
                
                Means_3d = zeros(1,1,size(Means.forecast,2));
                Means_3d(1,1,:) = Means.forecast(ind_r,:,peachcount+1);
                
                if fancharts
                    for i = 1:5
                        Bands.forecast(i,ind_r,:,peachcount+1) = min(Bands.forecast(i,ind_r,:,peachcount+1),Means_3d(1,1,:));
                    end
                    for i = 6:10
                        Bands.forecast(i,ind_r,:,peachcount+1) = max(Bands.forecast(i,ind_r,:,peachcount+1),Means_3d(1,1,:));
                    end
                else
                    Bands.forecast(1,ind_r,:,peachcount+1) = min(Bands.forecast(1,ind_r,:,peachcount+1),Means_3d(1,1,:));
                    Bands.forecast(2,ind_r,:,peachcount+1) = max(Bands.forecast(2,ind_r,:,peachcount+1),Means_3d(1,1,:));
                end
            end
        end
        
        Means.forecast = Means.forecast(1:ndim,:,:);
        Bands.forecast = Bands.forecast(:,1:ndim,:,:);
        
    end
    
    %% Four Quarter Forecasts and Q4Q4 Table
    
    % Calculate means and bands
    if any(ismember(plotList,{'4Q Forecast','Q4Q4 Table'}))
        for peachcount = 0:peachflag
            for Ivar = 1:ndim
                switch peachcount
                    case 1
                        yp_4q = getyp_4q(ypath_cond,qahead,Ivar,cum_for_4q,popadj,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                    case 0
                        yp_4q = getyp_4q(ypath,qahead,Ivar,cum_for_4q,popadj,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                end
                
                Means.forecast4q(Ivar,:,peachcount+1) = mean(yp_4q,1);
                if fancharts
                    Bands.forecast4q([1 10],Ivar,:,peachcount+1) = hpdint(yp_4q,0.90,0);
                    Bands.forecast4q([2 9],Ivar,:,peachcount+1) = hpdint(yp_4q,0.80,0);
                    Bands.forecast4q([3 8],Ivar,:,peachcount+1) = hpdint(yp_4q,0.70,0);
                    Bands.forecast4q([4 7],Ivar,:,peachcount+1) = hpdint(yp_4q,0.60,0);
                    Bands.forecast4q([5 6],Ivar,:,peachcount+1) = hpdint(yp_4q,0.50,0);
                else
                    Bands.forecast4q(:,Ivar,:,peachcount+1) = hpdint(yp_4q,percent,0);
                end
            end
            
            if zerobound
                for Ivar = 1:ndim
                    switch peachcount
                        case 1
                            yp_4q_bdd = getyp_4q(ypath_cond_bdd,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                        case 0
                            yp_4q_bdd = getyp_4q(ypath_bdd,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                    end
                    if fancharts
                        Bands.forecast4q([1 10],Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,0.90,0);
                        Bands.forecast4q([2 9],Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,0.80,0);
                        Bands.forecast4q([3 8],Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,0.70,0);
                        Bands.forecast4q([4 7],Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,0.60,0);
                        Bands.forecast4q([5 6],Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,0.50,0);
                    else
                        Bands.forecast4q(:,Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,percent,0);
                    end
                end
                
                ind_r = find(strcmp(varnames,'Interest Rate'));
                
                Means_3d_4q = zeros(1,1,size(Means.forecast4q,2));
                Means_3d_4q(1,1,:) = Means.forecast4q(ind_r,:,peachcount+1);
                
                if fancharts
                    for i = 1:5
                        Bands.forecast4q(i,ind_r,:,peachcount+1) = min(Bands.forecast4q(i,ind_r,:,peachcount+1),Means_3d_4q(1,1,:));
                    end
                    for i = 6:10
                        Bands.forecast4q(i,ind_r,:,peachcount+1) = max(Bands.forecast4q(i,ind_r,:,peachcount+1),Means_3d_4q(1,1,:));
                    end
                else
                    Bands.forecast4q(1,ind_r,:,peachcount+1) = min(Bands.forecast4q(1,ind_r,:,peachcount+1),Means_3d_4q(1,1,:));
                    Bands.forecast4q(2,ind_r,:,peachcount+1) = max(Bands.forecast4q(2,ind_r,:,peachcount+1),Means_3d_4q(1,1,:));
                end
                
            end
            
            if any(ismember(plotList,'Q4Q4 Table'))
                for Ivar = 1:ndim
                    switch peachcount
                        case 1
                            yp_4q = getyp_4q(ypath_cond,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                        case 0
                            yp_4q = getyp_4q(ypath,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                    end
                    
                    Means.table4q(Ivar,:,peachcount+1) = mean(yp_4q,1);
                    Bands.table4q(:,Ivar,:,peachcount+1) = hpdint(yp_4q,percent,0);
                    
                    if zerobound
                        switch peachcount
                            case 1
                                yp_4q_bdd = getyp_4q(ypath_cond_bdd,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                            case 0
                                yp_4q_bdd = getyp_4q(ypath_bdd,qahead,Ivar,cum_for_4q,popadj_plot,dlpop_fcst,YY,Idate,dlpop_old,q_adj);
                        end
                        
                        Bands.table4q(:,Ivar,:,peachcount+1) = hpdint(yp_4q_bdd,percent,0);
                        % If we report the interest rate forecast in Q4Q4 tables, should we use the means and bands as adjusted above or unadjusted?
                    end
                end
            end
        end
    end
    
    %% Counterfactual
    if any(ismember(plotList,{'Counterfactual by Variable','Counterfactual by Shock'}))
        
        for peachcount = 0:peachflag
            shockcount = 0;
            for Shockremove = ShockremoveList
                
                shockcount = shockcount + 1;
                
                % Read in data
                eval(['fid_ctr = fid1' num2str(peachcount) num2str(Shockremove) ';']);
                ypath_counter = loadcounter(fid_ctr,num_counter,ndim,counter_ahead,nsim,jstep);
                
                
                eval(['fclose(fid1' num2str(peachcount) num2str(Shockremove) ');']);
                
                % Calculate means and bands
                for Ivar = 1:ndim
                    yp_counter = getyp_counter(ypath_counter,counter_ahead,Ivar,cum_for_plot,popadj_plot,mnobss,dlpop_old,dlpop_ctr,repmat(YY(Startdate,:),size(ypath_counter,1),1),q_adj);
                    Means.counter(Ivar,:,peachcount+1,shockcount) = mean(yp_counter,1);
                    Bands.counter(:,Ivar,:,peachcount+1,shockcount) = hpdint(yp_counter,percent,0);
                end
                
                disp(['Loaded Counterfactual ' num2str(Shockremove)])
                
            end
            
            switch peachcount
                case 0, disp('Loaded All Counterfactuals (Unconditional)')
                case 1, disp('Loaded All Counterfactuals (Conditional)')
            end
            
        end
        clear ypath_counter yp_counter
    end
    
    %% Shock Decomposition
    if any(ismember(plotList,{'Shock Decomposition','Shock Decomposition-noDet'}))
        
        % Read in data
            ytrend = [];
            while(ftell(fid200) < num_trend)
                ypathadd = fread(fid200,[ndim,nsim/jstep],'single')';
                ytrend = [ytrend;ypathadd];
                ftell(fid200);
            end
            clear ypathadd
            fclose(fid200);
        
        
        ypath_trend = zeros(size(ytrend,1),counter_ahead*ndim);
        for Ivar = 1:ndim
            ypath_trend(:,counter_ahead*(Ivar-1)+1:counter_ahead*Ivar) = repmat(ytrend(:,Ivar),1,counter_ahead);
        end
        
        % Calculate means and bands
        Means_trend_raw = mean(ypath_trend,1);
        
        for Ivar = 1:ndim
            yp_trend = getyp_counter(ypath_trend,counter_ahead,Ivar,cum_for_plot,popadj_plot,mnobss,dlpop_old,dlpop_ctr,ytrend,q_adj);
            if exist('medianFlag','var') && medianFlag==1
                Means.trend(Ivar,:) = median(yp_trend,1);
            else
                Means.trend(Ivar,:) = mean(yp_trend,1);
            end
            
            Bands.trend(:,Ivar,:) = hpdint(yp_trend,percent,0);
        end
        
        clear ypath_trend yp_trend ytrend
        
        % The lines below load the deterministic trend.
        
        ypath_dettrend_all = loadcounter(fid300,num_shockdec,ndim,counter_ahead+1,nsim,jstep);
        % Explanation of the lines below: To turn log levels of
        % hours into a growth rate you need to use the previous
        % period's data. However, since the entire counterfactual
        % path (from t = 1 on) is different from the data this will
        % give you weird numbers. Instead you need to save and load
        % the previous period's counterfactual observables and use
        % those. This requires the manipulations below...
        
        yend = zeros(size(ypath_dettrend_all,1),ndim);
        ypath_dettrend = zeros(size(ypath_dettrend_all,1),ndim*counter_ahead);
        
        for Ivar = 1:ndim
            yend(:,Ivar) = ypath_dettrend_all(:,(Ivar-1)*(counter_ahead+1)+1);
            ypath_dettrend(:,(Ivar-1)*counter_ahead+1:Ivar*counter_ahead) = ypath_dettrend_all(:,(Ivar-1)*(counter_ahead+1)+2:Ivar*(counter_ahead+1));
        end
        
        fclose(fid300);
        
        if counter_ahead - mnobss > size(dlMA_pop,1)
            dlpop_ctr = [repmat(dlMA_pop',size(yend,1),1), repmat(dlMA_pop(end),size(yend,1),counter_ahead-mnobss-size(dlMA_pop,1))];
        else
            dlpop_ctr = repmat(dlMA_pop(1:counter_ahead-mnobss)',size(yend,1),1);
        end
        dlpop_old = repmat(dlpop(end-mnobss+1:end)',(size(yend,1)),1);
        
        for Ivar = 1:ndim
            yp_dettrend = getyp_counter(ypath_dettrend,counter_ahead,Ivar,cum_for_plot,popadj_plot,mnobss,dlpop_old,dlpop_ctr,yend,q_adj);
            Means.dettrend(Ivar,:,1) = mean(yp_dettrend,1);
        end
        
        % The lines below do the same for the determinisic trend smoothed
        % with peachdata.
        
        Means_shockdec_raw = zeros(length(ShockremoveList+1),size(ypath_dettrend,2));
        Means_shockdec_raw(end,:) = mean(ypath_dettrend,1);
        
        if peachflag
            ypath_dettrend_all = loadcounter(fid301,num_shockdec,ndim,counter_ahead+1,nsim,jstep);
            yend = zeros(size(ypath_dettrend_all,1),ndim);
            ypath_dettrend = zeros(size(ypath_dettrend_all,1),ndim*counter_ahead);
            
            for Ivar = 1:ndim
                yend(:,Ivar) = ypath_dettrend_all(:,(Ivar-1)*(counter_ahead+1)+1);
                ypath_dettrend(:,(Ivar-1)*counter_ahead+1:Ivar*counter_ahead) = ypath_dettrend_all(:,(Ivar-1)*(counter_ahead+1)+2:Ivar*(counter_ahead+1));
            end
            
            if counter_ahead - mnobss > size(dlMA_pop,1)
                dlpop_ctr = [repmat(dlMA_pop',size(yend,1),1), repmat(dlMA_pop(end),size(yend,1),counter_ahead-mnobss-size(dlMA_pop,1))];
            else
                dlpop_ctr = repmat(dlMA_pop(1:counter_ahead-mnobss)',size(yend,1),1);
            end
            dlpop_old = repmat(dlpop(end-mnobss+1:end)',(size(yend,1)),1);
            
            for Ivar = 1:ndim
                yp_dettrend = getyp_counter(ypath_dettrend,counter_ahead,Ivar,cum_for_plot,popadj_plot,mnobss,dlpop_old,dlpop_ctr,yend,q_adj);
                Means.dettrend(Ivar,:,2) = mean(yp_dettrend,1);
            end
        end
        
        clear ypath_dettrend_all ypath_dettrend yp_dettrend
        
        for peachcount = 0:peachflag
            
            shockcount = 0;
            
            for Shockremove = ShockremoveList
                
                shockcount = shockcount + 1;
                
                eval(['fid_ctr = fid2' num2str(peachcount) num2str(Shockremove) ';']);
                ypath_shockdec_all = loadcounter(fid_ctr,num_shockdec,ndim,counter_ahead+1,nsim,jstep);
                yend = zeros(size(ypath_shockdec_all,1),ndim);
                ypath_shockdec = zeros(size(ypath_shockdec_all,1),ndim*counter_ahead);
                
                for Ivar = 1:ndim
                    yend(:,Ivar) = ypath_shockdec_all(:,(Ivar-1)*(counter_ahead+1)+1);
                    ypath_shockdec(:,(Ivar-1)*counter_ahead+1:Ivar*counter_ahead) = ypath_shockdec_all(:,(Ivar-1)*(counter_ahead+1)+2:Ivar*(counter_ahead+1));
                end
                
                Means_shockdec_raw(shockcount,:) = mean(ypath_shockdec,1);
                
                eval(['fclose(fid2' num2str(peachcount) num2str(Shockremove) ');']);
                
                if counter_ahead - mnobss > size(dlMA_pop,1)
                    dlpop_ctr = [repmat(dlMA_pop',size(yend,1),1), repmat(dlMA_pop(end),size(yend,1),counter_ahead-mnobss-size(dlMA_pop,1))];
                else
                    dlpop_ctr = repmat(dlMA_pop(1:counter_ahead-mnobss)',size(yend,1),1);
                end
                dlpop_old = repmat(dlpop(end-mnobss+1:end)',(size(yend,1)),1);
                
                for Ivar = 1:ndim
                    yp_shockdec = getyp_counter(ypath_shockdec,counter_ahead,Ivar,cum_for_plot,popadj_plot,mnobss,dlpop_old,dlpop_ctr,yend,q_adj);
                    
                    if exist('medianFlag','var') && medianFlag==1
                        Means.shockdec(Ivar,:,peachcount+1,shockcount) = median(yp_shockdec,1);
                    else
                        Means.shockdec(Ivar,:,peachcount+1,shockcount) = mean(yp_shockdec,1);
                    end
                    
                    Bands.shockdec(:,Ivar,:,peachcount+1,shockcount) = hpdint(yp_shockdec,percent,0);
                end
                
                disp(['Loaded Shock Decomposition ' num2str(Shockremove)])
                
            end
            
            switch peachcount
                case 0, disp('Loaded All Shock Decompositions (Unconditional)')
                case 1, disp('Loaded All Shock Decompositions (Conditional)')
            end
            
        end
        
        clear ypath_shockdec yp_shockdec ypath_shockdec_all yend
        
    end
    
    %% Shock
    
    if any(ismember(plotList,'Shock'))
        
        % Read in data
        datashocks = [];
        dickdatashocks = [];
        semidickdatashocks = [];
        while ( ftell(fid4) < num_data )
            datashocksadd = fread(fid4,[nshocks*nobs,nsim/jstep],'single')';
            datashocks = [datashocks;datashocksadd];
            
            ftell(fid4);
            
            if peachflag
                dickdatashocksadd = fread(fid3,[nshocks*nobs,nsim/jstep],'single')';
                dickdatashocks = [dickdatashocks;dickdatashocksadd];
                clear dickdatashocksadd;
                ftell(fid3);
            end
        end
        if peachflag
            fclose(fid3);
        end
        fclose(fid4);
        
        if peachflag
            dickspath = [];
            while ( ftell(fid2) < num_forshocks )
                dickspathadd = fread(fid2,[nshocks*psize,nsim/jstep],'single')';
                dickspath = [dickspath;dickspathadd];
                ftell(fid2);
            end;
            fclose(fid2);
            
        end
        
        % Calculate means and bands
        for peachcount = 0:peachflag
            for III = 1:nshocks
                switch peachcount
                    case 1
                        dshockfor = dickspath(:,psize*(III-1)+1:psize*III);
                        means_dickshocks(III,:) = mean(dshockfor,1);
                        bands_dickshocks(:,III,:,:) = hpdint(dshockfor,percent,0);
                        
                        dshocksample = dickdatashocks(:,III:nshocks:end-mod(size(dickdatashocks,2),nshocks));
                        means_dickdatashocks(III,:) = mean(dshocksample,1);
                        bands_dickdatashocks(:,III,:,:) = hpdint(dshocksample,percent,0);
                    case 0
                        shocksample = datashocks(:,III:nshocks:end-mod(size(datashocks,2),nshocks));
                        means_datashocks(III,:) = mean(shocksample,1);
                        bands_datashocks(:,III,:,:) = hpdint(shocksample,percent,0);
                end
            end
            
            switch peachcount,
                case 1
                    Means.condshocks = [means_dickdatashocks means_dickshocks];
                    Bands.condshocks = cat(3,bands_dickdatashocks,bands_dickshocks);
                case 0
                    Means.uncondshocks = means_datashocks;
                    Bands.uncondshocks = bands_datashocks;
            end
        end
        
    end
    if any(ismember(plotList,'ShockAnt'))
        
        % Read in data
        datashocks_ns = [];
        dickdatashocks_ns = [];
        semidickdatashocks_ns = [];
        while ( ftell(fid45) < num_data )
            datashocksadd_ns = fread(fid45,[nshocks*nobs,nsim/jstep],'single')';
            datashocks_ns = [datashocks_ns;datashocksadd_ns];
            
            ftell(fid45);
            
            if peachflag
                dickdatashocksadd_ns = fread(fid35,[nshocks*nobs,nsim/jstep],'single')';
                dickdatashocks_ns = [dickdatashocks_ns;dickdatashocksadd_ns];
                clear dickdatashocksadd_ns;
                ftell(fid35);
            end
        end
        if peachflag
            fclose(fid35);
            
        end
        fclose(fid45);
        
        if peachflag
            dickspath_ns = [];
            while ( ftell(fid25) < num_forshocks )
                dickspathadd_ns = fread(fid25,[nshocks*psize,nsim/jstep],'single')';
                dickspath_ns = [dickspath_ns;dickspathadd_ns];
                ftell(fid25);
            end;
            fclose(fid25);
            
        end
        
        % Calculate means and bands
        for peachcount = 0:peachflag
            for III = 1:nshocks
                switch peachcount
                    case 1
                        dshockfor_ns = dickspath_ns(:,psize*(III-1)+1:psize*III);
                        means_dickshocks_ns(III,:) = mean(dshockfor_ns,1);
                        bands_dickshocks_ns(:,III,:,:) = hpdint(dshockfor_ns,percent,0);
                        
                        dshocksample_ns = dickdatashocks_ns(:,III:nshocks:end-mod(size(dickdatashocks_ns,2),nshocks));
                        means_dickdatashocks_ns(III,:) = mean(dshocksample_ns,1);
                        bands_dickdatashocks_ns(:,III,:,:) = hpdint(dshocksample_ns,percent,0);
                    case 0
                        shocksample_ns = datashocks_ns(:,III:nshocks:end-mod(size(datashocks_ns,2),nshocks));
                        means_datashocks_ns(III,:) = mean(shocksample_ns,1);
                        bands_datashocks_ns(:,III,:,:) = hpdint(shocksample_ns,percent,0);
                end
            end
            
            switch peachcount,
                case 1
                    Means.condshocks_ns = [means_dickdatashocks_ns means_dickshocks_ns];
                    Bands.condshocks_ns = cat(3,bands_dickdatashocks_ns,bands_dickshocks_ns);
                case 0
                    Means.uncondshocks_ns = means_datashocks_ns;
                    Bands.uncondshocks_ns = bands_datashocks_ns;
            end
        end
        
    end
    
    % Save means and bands
    if overwrite, save([spath,'/fcastMeans'],'Means','Bands'); end
    
else
    
    load([spath,'/fcastMeans.mat'],'Means','Bands');
end

fprintf('Means and bands are saved in %s. \n', spath)


