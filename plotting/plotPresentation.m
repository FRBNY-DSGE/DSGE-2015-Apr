% OVERVIEW (WRITTEN AND INTRODUCED ON FEB 18, 2014)
% plotPresentation.m: Produces figures for presentations only.
%
% BIG UPDATE: Can plot either states or observables depending on the
% value of vars_to_plot, which is set in forplot.m and forplot_mode.m.
% This replaces plotPresentation_states.m and the old plotPresentation.m
%
% IMPORTANT FUNCTIONS
% setVseq.m: a function whose input is the figure type. The output from
% this
%            function is a vector that indicates the sequence of variables
%            for which we want to produce figures.
% recordFcn.m: a function that records the plotted values.
%              The output is a structure variable, record, whose fields
%              are different plot types, and whose subfields are the
%              plotted
%              mean and bands. The variable record is saved and
%              written to an output text file.
% figspecs.m: this function outputs the variables Xaxis, Yaxis, Title, line,
%             lgnd, and plotSeperate, which are all specifications
%             for unique to each figure and model. lgnd contains
%             specifications for the legend. plotSeparate is a flag
%             indicating whether or not the figures for each variable are
%             to be plotted separately.
% setfig.m: this function takes the variable plotType (ie. 'Forecast') as
%           input and outputs the scalar variable figmain. All figures of
%           this plot type have a figure number that is a function of this
%           scalar. For example, if figmain for 'Forecast' figures is 100,
%           then the output forecast figure would be figure(101), and
%           figure(102) for labor share forecast.
% applyfigspecs.m: this program applies figure specifications for the
%                  X-axis, Y-axis, title, secondary axis, and gridlines.

% Set up 'number of' dimension to either nvar or nplotstates
    ndim = nvar;
    fcast_suffix = '';
    cum_for_plot = cum_for;
    popadj_plot = popadj;


pres = 1; newsletter = 0; system = 0;
peachloop = [1];
noTitles = 0; % 1 for no titles on the plots


%% Specify outfile for recording plotted values
outfile0 = [gpath,date,'.txt'];
fid0 = fopen(outfile0,'wt+');

%% Transform actual data

    Ys = getYs(YY(1:Idate,:),YY_p,dlpop,dlpop_p,cum_for_plot,popadj_plot,ndim,fourq_flag,q_adj,0);


    if any(ismember(plotList,{'Counterfactual by Variable','Counterfactual by Shock','Shock Decomposition'}))
        Ys_counter = Ys(1:Startdate,:);
    end


%% Loop through figures
for plotnum = 1:length(plotList)
    plotType = plotList{plotnum};
    
    sirf = Startdate:stime+future;
    
    [Vseq,Vseq_alt,Shseq,Shseq_alt] = setVseq(plotType,nshocks,nant,ndim,mspec);


    for peachcount = peachloop
        varcount = 0;
        V_alt_count = 0;
        for V_idx = 1:length(Vseq)

            V_a = Vseq(V_idx);
            pass = 0;

            %% Set figure specifications and fignum
            [Xaxis,Yaxis,Title,line,lgnd,plotSeparate] = feval(['figspecs', fcast_suffix], ...
                V_a,V_idx,mspec,peachcount,plotType,...
                varnames,varnames_YL,varnames_YL_4Q,Vseq,Vseq_alt,nobs,Idate,Startdate,Enddate,...
                sirf,sirf_shockdec,sirf_counter,sirf_shock_1,datesall,dataset,zerobound,...
                ndim,list,Counterforecast,shocksnames,names_shocks_title,pass);
            
            % This extends the bounds for graphs with the (larger) 99
            % percent bands
            if onepctflag && strcmp(plotType,'Forecast')
                if any(V_a == [1,2])
                    Yaxis.lim = [-15,15];
                elseif any(V_a == [4,5])
                    Yaxis.lim = [-5,10];
                end
            end
            
            [figmain] = setfig(plotType);
            
            %% Specify subplots and orientation
            
            if plotSeparate
                figure(figmain+(250*peachcount)+V_a);
            else
                if strcmp(plotType,'Shock')
                    figure(figmain+(250*peachcount))
                end
                switch plotType
                    case {'Shock Decomposition','Forecast','Counterfactual by Variable','Counterfactual by Shock'} % Figures that loop through variables
                        if any(strcmp(plotType,{'Forecast','Counterfactual by Shock'}))
                            figure(figmain); % Both conditional and unconditional forecasts are plotted together
                            %subplot(length(Vseq)-length(Vseq_alt),2,2*V_idx-(1-round(peachcount/2)));
                        elseif strcmp(plotType,'Shock Decomposition')
                            figure(figmain+(250*peachcount));
                            orient portrait
                            ax = get(gca,'Position');
                            set(gca,'Position',ax);
                        end
                end
            end
            
            %% Plot figures
            switch plotType
                case {'Forecast'}
                    % Mean
                    mlp = [Ys(:,V_a);squeeze(Means.forecast(V_a,:,peachcount+1))'];
                    
                    % Bands
                    num_bands = 2 + 8*fancharts;
                    bandsp = [repmat(Ys(:,V_a),[1 num_bands]);squeeze(Bands.forecast(:,V_a,:,peachcount+1))'];
                    
                    LBands(V_a,:) =  bandsp(sirf,1)';
                    UBands(V_a,:) =  bandsp(sirf,2)';
                    
                    % Record data plotted
                    if overwrite
                        [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                            bandsp,mlp);
                    end
                    
                    % Fancharts
                    if fancharts
                        bandsp = [bandsp(:,1),bandsp(:,2:10)-bandsp(:,1:9)];
                    end
                    
                    if fancharts
                        if mspec==990 && isnan(bandsp(stime,1))
                                P = area(stime+1:sirf(end),bandsp(stime+1:sirf(end),:));      
                        else
                            P = area(sirf,bandsp(sirf,:));
                        end

                        set(P(1),'FaceColor','none')
                        for layer = 2:10
                            set(P(layer),'FaceColor',[0.4 + 0.12*abs(layer-6),0.4 + 0.12*abs(layer-6),1]);
                        end
                        set (P,'LineStyle','none')
                        if V_a == 3
                            set(P,'BaseValue',.52)
                        end
                    elseif ~strcmp(parstr,'Mode')
                        P = plot(sirf,bandsp(sirf,1),'r:',sirf,bandsp(sirf,2),'r:');
                        set(P,'LineWidth',line.width);
                    end
                    hold on;
                    
                    % This adds the 99 percent bands if onepctflag is on
                    if onepctflag
                        bandsp_onepct = [repmat(Ys(:,V_a),1,2);squeeze(Bands_onepct(:,V_a,:,peachcount+1))'];
                        P = plot(sirf,bandsp_onepct(sirf,1),'r:',sirf,bandsp_onepct(sirf,2),'r:');
                        set(P,'LineWidth',line.width);
                    end
                    
                    % Forecast
                    P = plot(sirf,mlp(sirf),'r-');
                    set(P,'LineWidth',line.width);
                    
                    hold on;
                    
                    % Actual data
                    mlp = [Ys(:,V_a);NaN*squeeze(Means.forecast(V_a,:,peachcount+1))'];
                    P = plot(sirf,mlp(sirf),'k-');
                    set(P,'LineWidth',line.width);
                    
                    hold off;
                    
                case 'Shock Decomposition'
                    
                    % Means_cat
                    eval(strcat('states',num2str(mspec)));
                    Means_cat = zeros(ndim,counter_ahead,length(shockcats));
                    for j = 1:length(shockcats)
                        catmean = zeros(size(Means.trend));
                        catshocks = shockcats{j};
                        for i2 = 1:length(shockcats{j})
                            % Use 0 to indicate deterministic trend
                            if catshocks(i2) == 0
                                catmean = catmean + Means.dettrend(:,:,peachcount+1) - Means.trend;  %--> we get Det.Trend for each observable (for each period)
                                % All other shocks...
                            else
                                catmean = catmean + Means.shockdec(:,:,peachcount+1,catshocks(i2)) - Means.trend;  %--> we get the other "categorical" shocks
                            end
                        end
                        Means_cat(:,:,j) = catmean(1:ndim,:);
                    end
                    if any(strcmp('Residual',list))
                        ind_res = find(strcmp('Residual',list));
                        Means_cat(:,:,ind_res) = Means_cat(:,:,ind_res) + Means.dettrend(:,:,peachcount+1) - Means.trend;
                    end
                   
                    % Shock decompositions
                    if size(shockcats,1)==1
                        mlp_bar = [NaN(Startdate,size(Means_cat,3)); squeeze(Means_cat(V_a,:,:))'];
                    else
                        mlp_bar = [NaN(Startdate,size(Means_cat,3)); squeeze(Means_cat(V_a,:,:))];
                    end
                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % ---> mlp_bar gives us values for shock decomposition bars
                    barh = bar(sirf_shockdec,mlp_bar(sirf_shockdec,:),1.5);
                    % set bar colors
                    for iList = 1:length(list)
                        if isnumeric(shockdec_color{iList})
                            set(barh(iList),'FaceColor',shockdec_color{iList});
                        else
                            set(barh(iList),'FaceColor',rgb(shockdec_color{iList}));
                        end
                    end
                    
                    set(barh,'EdgeColor','none')
                    hold on;
                                        
                    % Forecast
                    mlp = [Ys(:,V_a);squeeze(Means.forecast(V_a,:,peachcount+1))'];
                    mlp = mlp - [NaN(1,size(mlp,1) - size(Means.trend,2)) squeeze(Means.trend(V_a,:))]';
                    %---> mlp gives us detrended data + forecast (i.e. we get the deviations from steady state here)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Record data plotted
                    if overwrite
                        [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                            list,mlp,mlp_bar);
                    end
                    
                    P = plot(sirf_shockdec,mlp(sirf_shockdec),'r-');
                    set(P,'LineWidth',line.width);
                    
                    % Actual data
                    mlp = [Ys(:,V_a);NaN*squeeze(Means.forecast(V_a,:,peachcount+1))'];
                    mlp = mlp - [NaN(1,size(mlp,1) - size(Means.trend,2)) squeeze(Means.trend(V_a,:))]'; % This line detrends the data/forecast
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    P = plot(sirf_shockdec,mlp(sirf_shockdec),'k-');
                    set(P,'LineWidth',line.width);
                    
                    % Legend
                    legh = legend(list,'Location','SouthOutside','Orientation','horizontal');


                    legpos = get(legh,'Position');
                    plotpos = get(gca,'Position');
                    legpos(1) = plotpos(1) + plotpos(3)/2 - legpos(3)/2;
                    legpos(2) = plotpos(2) - plotpos(4)/6 - legpos(4);
                        
                    set(legh,'Position',legpos);
                    
                case 'Counterfactual by Variable'
                    
                    
                    for S_1 = 1:length(Shseq)
                        S_a = Shseq(S_1);
                        % to specify title for graphs
                        pass=1;
                        [Xaxis,Yaxis,Title,line,lgnd,plotSeparate] = ...
                          figspecs(V_a,V_idx,mspec,peachcount,plotType,...
                            varnames,varnames_YL,varnames_YL_4Q,Vseq,Vseq_alt,nobs,Idate,Startdate,Enddate,...
                            sirf,sirf_shockdec,sirf_counter,sirf_shock_1,datesall,dataset,zerobound,...
                            ndim,list,Counterforecast,shocksnames,names_shocks_title,pass,S_a);
                        
                        if plotSeparate
                            figure(figmain+(250*peachcount)+(10*V_a)+S_a);
                        else
                            thisfignum=figmain+(250*peachcount)+(10*V_a);
                            if S_1>length(Shseq)/2;
                                thisfignum=thisfignum+1;
                            end
                            figure(thisfignum);
                            %set(gcf, 'PaperSize', [10 40]);
                            rows=ceil(length(Shseq))/4;
                            cols=2;
                            item_num=mod(S_1,length(Shseq)/2);
                            if item_num==0; item_num=length(Shseq)/2; end;
                            s=subplot(rows,cols, item_num);
                            if rem(S_1,2)==1 % odd
                                col=1;
                            else             % even
                                col=2;
                            end
                            row=ceil(S_1/2);
                            %axes('Position',[(col-1)/cols 1-(row-1)/rows 1/cols 1/rows]);
                        end
                        
                        % Counterfactual forecast
                        mlp_counter_byvar = [Ys_counter(:,V_a);squeeze(Means.counter(V_a,:,peachcount+1,S_a))'];
                        bandsp_counter_byvar = [repmat(Ys_counter(:,V_a),[1 2]);squeeze(Bands.counter(:,V_a,:,peachcount+1,S_a))'];
                        LBands_counter_byvar(V_a,:) =  bandsp_counter_byvar(sirf_counter,1)';
                        UBands_counter_byvar(V_a,:) =  bandsp_counter_byvar(sirf_counter,2)';
                        
                        P = plot(sirf_counter,mlp_counter_byvar(sirf_counter),'g-');
                        set(P,'LineWidth',line.width);
                        hold on;
                        
                        % Record data plotted
                        if overwrite
                            [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                                S_1,S_a,Shseq,mlp_counter_byvar);
                        end
                        orient tall
                        set(gcf, 'PaperPosition', [0 0 8 10.5]);
                                                
                        % Forecast
                        mlp_counter_byvar = [Ys_counter(:,V_a);NaN*squeeze(Means.counter(V_a,:,peachcount+1,S_a))'];
                        mlp = [Ys(:,V_a);squeeze(Means.forecast(V_a,:,peachcount+1))'];
                        
                        P = plot(sirf_counter,mlp(sirf_counter),'r-');
                        set(P,'LineWidth',line.width);
                        
                        % Actual Data
                        mlp = [Ys(:,V_a);NaN*squeeze(Means.forecast(V_a,:,peachcount+1))'];
                        P = plot(sirf_counter,mlp(sirf_counter),'k-');
                        set(P,'LineWidth',line.width);
                        
                        % Figure Specifications
                        applyfigspecs;
                    end
                    
                case 'Counterfactual by Shock'
                    
                    % NOTE: the first loop is Vseq, but these are actually the shocks that we want to remove.
                    % The second loop is Shseq, but these are the variables that we are plotting counterfactuals for.
                    
                    shockcount = 0; % This actually counts the number of variables (not shocks) we've looped through
                    for S_1 = 1:length(Shseq)
                        S_a = Shseq(S_1);
                        
                        pass=1;
                        [Xaxis,Yaxis,Title,line,lgnd,plotSeparate] = ...
                          figspecs(V_a,V_idx,mspec,peachcount,plotType,...
                            varnames,varnames_YL,varnames_YL_4Q,Vseq,Vseq_alt,nobs,Idate,Startdate,Enddate,...
                            sirf,sirf_shockdec,sirf_counter,sirf_shock_1,datesall,dataset,zerobound,...
                            ndim,list,Counterforecast,shocksnames,names_shocks_title,pass,S_a);
                        
                        if plotSeparate_counter
                            fignum = figmain + 250*peachcount + 10*S_a + V_a;
                            figure(fignum);
                        else
                            figure(figmain+250*peachcount+10*V_a);
                            %set(gcf, 'PaperSize', [10 40]);
                            s=subplot(ceil(length(Shseq)/2),2,S_1);
                        end
                        
                        
                        % Counterfactual forecast
                        mlp_counter = [Ys_counter(:,S_a);squeeze(Means.counter(S_a,:,peachcount+1,V_a))'];
                        bandsp_counter = [repmat(Ys_counter(:,S_a),[1 2]);squeeze(Bands.counter(:,S_a,:,peachcount+1,V_a))'];
                        LBands_counter(S_a,:) =  bandsp_counter(sirf_counter,1)';
                        UBands_counter(S_a,:) =  bandsp_counter(sirf_counter,2)';
                        
                        P = plot(sirf_counter,mlp_counter(sirf_counter),'g-');
                        set(P,'LineWidth',line.width);
                        hold on;
                        
                        % Record data plotted
                        if overwrite
                            [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                                S_1,S_a,Shseq,mlp_counter);
                        end

                        orient tall
                        % Forecast
                        mlp_counter = [Ys_counter(:,S_a);NaN*squeeze(Means.counter(S_a,:,peachcount+1,V_a))'];
                        mlp = [Ys(:,S_a);squeeze(Means.forecast(S_a,:,peachcount+1))'];
                        
                        P = plot(sirf_counter,mlp(sirf_counter),'r-');
                        set(P,'LineWidth',line.width);
                        
                        % Actual data
                        mlp = [Ys(:,S_a);NaN*squeeze(Means.forecast(S_a,:,peachcount+1))'];
                        P = plot(sirf_counter,mlp(sirf_counter),'k-');
                        set(P,'LineWidth',line.width);
                        
                        % Figure Specifications
                        applyfigspecs;
                    end
                    
                case 'Shock'
                    sirf_shock = [sirf_shock_1,(repmat(Idate,1,pnobss))+(logical(peachcount)*(1:pnobss))]; % if conditional, adjust sirf_shock to pnobss (psize) periods
                    
                    % Dotted line at zero
                    P = plot(sirf_shock,zeros(size(sirf_shock)),'k--');
                    hold on;
                    
                    % Shock path
                    if peachcount
                        if peachcount == 2
                            mlp = squeeze(Means.semicondshocks(V_a,:))';
                            bandsp = squeeze(Bands.semicondshocks(:,V_a,:))';
                        else
                            mlp = squeeze(Means.condshocks(V_a,:))';
                            bandsp = squeeze(Bands.condshocks(:,V_a,:))';
                        end
                    else
                        mlp = squeeze(Means.uncondshocks(V_a,:))';
                        bandsp = squeeze(Bands.uncondshocks(:,V_a,:))';
                    end
                    P = plot(sirf_shock,mlp(sirf_shock),'r-');
                    set(P,'LineWidth',line.width);
                    hold on;
                    
                    % Record data plotted
                    if overwrite
                        [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                            Vseq,mlp,pnobss);
                    end
                    
                    % In sample
                    insample = find(sirf_shock<nobs+1);
                    outsample = find(sirf_shock>nobs);
                    if ~isempty(insample)
                        if peachcount

                                mlp = [squeeze(Means.condshocks(V_a,sirf_shock(insample)))';NaN*squeeze(Means.condshocks(V_a,sirf_shock(outsample)))'];
                        else
                            mlp = [squeeze(Means.uncondshocks(V_a,sirf_shock(insample)))';NaN*squeeze(Means.uncondshocks(V_a,sirf_shock(outsample)))'];
                        end
                        P = plot(sirf_shock',mlp,'k-');
                        set(P,'LineWidth',line.width);
                    end;
                    
                    P = plot(sirf_shock,bandsp(sirf_shock,1),'k:',sirf_shock,bandsp(sirf_shock,2),'k:');

                case 'ShockAnt'
                    sirf_shock = [sirf_shock_1,(repmat(Idate,1,pnobss))+(logical(peachcount)*(1:pnobss))]; % if conditional, adjust sirf_shock to pnobss (psize) periods
                    
                    % Dotted line at zero
                    P = plot(sirf_shock,zeros(size(sirf_shock)),'k--');
                    hold on;
                    
                    % Shock path
                    switch peachcount
                        case 0, 
                            mlp = squeeze(Means.uncondshocks_ns(V_a,:))';
                            bandsp = squeeze(Bands.uncondshocks_ns(:,V_a,:))';
                        case 1, 
                            mlp = squeeze(Means.condshocks_ns(V_a,:))';
                            bandsp = squeeze(Bands.condshocks_ns(:,V_a,:))';
                        case 2, 
                            mlp = squeeze(Means.semicondshocks_ns(V_a,:))';
                            bandsp = squeeze(Bands.semicondshocks_ns(:,V_a,:))';
                    end
                    P = plot(sirf_shock,mlp(sirf_shock),'r-');
                    set(P,'LineWidth',line.width);
                    hold on;

                    % Record data plotted
                    if overwrite
                        [record] = recordFcn(plotType,record,V_idx,V_a,peachcount,peachfile,tiall,Startdate,stime,qahead,varnames,names_shocks,fid0,0,...
                                                Vseq,mlp,pnobss);
                    end

                    % In sample
                    insample = find(sirf_shock<nobs+1);
                    outsample = find(sirf_shock>nobs);
                    if ~isempty(insample)
                        switch peachcount
                            case 1, mlp = [squeeze(Means.condshocks_ns(V_a,sirf_shock(insample)))';NaN*squeeze(Means.condshocks_ns(V_a,sirf_shock(outsample)))'];
                            case 0, mlp = [squeeze(Means.uncondshocks_ns(V_a,sirf_shock(insample)))';NaN*squeeze(Means.uncondshocks_ns(V_a,sirf_shock(outsample)))'];
                            case 2, mlp = [squeeze(Means.semicondshocks_ns(V_a,sirf_shock(insample)))';NaN*squeeze(Means.semicondshocks_ns(V_a,sirf_shock(outsample)))'];
                        end
                        P = plot(sirf_shock',mlp,'k-');
                        set(P,'LineWidth',line.width);
                    end;

                    P = plot(sirf_shock,bandsp(sirf_shock,1),'k:',sirf_shock,bandsp(sirf_shock,2),'k:');
                    set(P,'LineWidth',line.width);
                    applyfigspecs;

            end
            
            if any(strcmp(plotType,{'Forecast','Shock Decomposition','Shock'}))
                applyfigspecs; % For the counterfactuals (which have another loop in addition to the variable loop),
                % applyfigspecs is called within the loop so
                % individual figures can be modified
            end
        end
    end
end

%% Record
% Save plotted values in .txt and .mat
if overwrite
    finalRecord = 1;
    [record] = recordFcn('',record,V_idx,V_a,peachcount,peachfile,...
        tiall,Startdate,stime,qahead,...
        varnames,names_shocks,fid0,finalRecord);
   
    %% Save figures
    % Removed lambdas and plotNoShocks from figure names
    
    for plotnum = 1:length(plotList)
        plotType = plotList{plotnum};
        if ~strcmp(plotType,'Q4Q4 Table')
            
            [figmain] = setfig(plotType);
                [Vseq,Vseq_alt,Shseq,Shseq_alt] = setVseq(plotType,nshocks,nant,ndim,mspec);

            
            switch plotType
                case 'Forecast'
                    for peachcount = peachloop
                        for V_a = Vseq
                            fignum = figmain+(250*peachcount)+V_a;
                            print(figure(fignum),'-depsc', [gpath,'forecast',peachstr{peachcount+1},char(graph_title(V_a)),'.eps']);
                            saveas(figure(fignum),[gpath,'forecast',peachstr{peachcount+1},char(graph_title(V_a))],'pdf');
                        end
                    end
                    
                case 'Shock Decomposition'
                    for peachcount = peachloop
                        if plotSeparate
                            for V_a = Vseq
                                fignum = figmain+(250*peachcount)+V_a;
                                print(figure(fignum),'-depsc', [gpath,'shockdec',peachstr{peachcount+1},char(graph_title(V_a)),'.eps']);
                                saveas(figure(fignum),[gpath,'shockdec',peachstr{peachcount+1},char(graph_title(V_a))],'pdf');
                            end
                        else
                        end
                    end
                    
                case {'Counterfactual by Shock','Counterfactual by Variable'}
                    %temp fix
                    plotSeparate = 0;
                    if plotSeparate
                        for peachcount = peachloop
                            for V_a = Vseq
                                for S_a = Shseq
                                    if strcmp(plotType,'Counterfactual by Shock')
                                        fignum = figmain+(250*peachcount)+(10*S_a)+V_a;
                                    elseif strcmp(plotType,'Counterfactual by Variable')
                                        fignum = figmain+(250*peachcount)+(10*V_a)+S_a;
                                    end
                                    print(figure(fignum),'-depsc',[gpath,'counter',peachstr{peachcount+1},char(names_shocks_title(S_a)),char(graph_title(V_a)),'.eps']);
                                    saveas(figure(fignum),[gpath,'counter',peachstr{peachcount+1},char(names_shocks_title(S_a)),char(graph_title(V_a))],'pdf');
                                end
                            end
                        end
                    else
                        if strcmp('Counterfactual by Shock',plotType)
                            for peachcount = peachloop
                                for S_a = Shseq
                                    fignum = figmain+(250*peachcount)+(10*S_a);
                                    print(figure(fignum),'-depsc',[gpath,'counter_shock',peachstr{peachcount+1},char(names_shocks_title(S_a)),'.eps']);
                                    saveas(figure(fignum),[gpath,'counter_shock',peachstr{peachcount+1},char(names_shocks_title(S_a))],'pdf');
                                end
                            end
                            
                        elseif strcmp('Counterfactual by Variable',plotType)
                            for peachcount = peachloop
                                for V_a = Vseq
                                    fignum = figmain+(250*peachcount)+(10*V_a);
                                    fignum2= fignum+1;
                                    print(figure(fignum),'-depsc',[gpath,'counter_var',peachstr{peachcount+1},char(graph_title(V_a)),'.eps']);
                                    print(figure(fignum2),'-depsc',[gpath,'counter_var',peachstr{peachcount+1},char(graph_title(V_a)),'_continued','.eps']);
                                    saveas(figure(fignum),[gpath,'counter_var',peachstr{peachcount+1},char(graph_title(V_a))],'pdf');
                                    saveas(figure(fignum2),[gpath,'counter_var',peachstr{peachcount+1},char(graph_title(V_a)) '_continued'],'pdf');
                                end
                            end
                        end
                    end
                case 'Shock'
                    for peachcount = peachloop
                        plotSeparate=1;
                        if plotSeparate
                            for V_a = Vseq
                                fignum = figmain+(250*peachcount)+V_a;
                                print(figure(fignum),'-depsc', [gpath,'shock',peachstr{peachcount+1},char(names_shocks_title(V_a)),'.eps']);
                                saveas(figure(fignum),[gpath,'shock',peachstr{peachcount+1},char(names_shocks_title(V_a))],'pdf');
                            end
                        else
                            fignum = figmain+(250*peachcount);
                            print(figure(fignum),'-depsc', [gpath,'shock',peachstr{peachcount+1},'.eps']);
                            saveas(figure(fignum),[gpath,'shock',peachstr{peachcount+1}],'pdf');
                        end
                    end

                case 'ShockAnt'
                    for peachcount = peachloop
                        plotSeparate = 1;
                        if plotSeparate
                            for V_a = Vseq
                                fignum = figmain+(250*peachcount)+V_a;
                                print(figure(fignum),'-depsc', [gpath,'shockAnt',peachstr{peachcount+1},'_',char(names_shocks_title(V_a)),'.eps']);
                                saveas(figure(fignum),[gpath,'shockAnt',peachstr{peachcount+1},char(names_shocks_title(V_a))],'pdf');
                            end
                        else
                            fignum = figmain+(250*peachcount);
                            print(figure(fignum),'-depsc', [gpath,'shockAnt',peachstr{peachcount+1},'_','.eps']); 
                            saveas(figure(fignum),[gpath,'shockAnt',peachstr{peachcount+1}],'pdf');
                        end
                    end
                    
            end
        end
    end
    fclose(fid0);

    curPath = pwd; 
    eval(['cd ',curPath]); 
end

fprintf('Figures are saved in gpath (%s) \n',gpath);
