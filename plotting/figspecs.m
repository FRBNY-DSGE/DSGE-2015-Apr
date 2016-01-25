% figspecs.m will output specifications for each type of product and figure

function [Xaxis,Yaxis,Title,line,lgnd,plotSeparate] = ...
  figspecs(V_a,V_1,mspec,peachcount,plotType,...
    varnames,varnames_YL,varnames_YL_4Q,Vseq,Vseq_alt,nobs,Idate,Startdate,Enddate,...
    sirf,sirf_shockdec,sirf_counter,sirf_shock,datesall,dataset,zerobound,...
    nvar,list,Counterforecast,shocksnames,names_shocks,pass,varargin)

if pass % For counterfactuals, pass=1 only when figspecs is called the second time; for counterfactuals, figspecs requires S_a as well as V_a (shock and var loop)
    S_a = varargin{1};
end

switch plotType
    %% Forecast
    case {'Forecast','Forecast Semicond', 'Forecast Comparison','4Q Forecast','Exp_Forecast','Forecast Comparison Separate','4Q Forecast Separate','Forward Guidance'}
        plotSeparate = 0;
        
        Xaxis.limits = [min(sirf),max(sirf)];
        
        if strcmp('4Q Forecast',plotType)
            Yaxis.axislabel = varnames_YL_4Q(V_a);
        else
            Yaxis.axislabel = varnames_YL(V_a);
        end
        
        if any(V_a==Vseq_alt)
            Title.name = [];
            Title.size = 13;
            Xaxis.freq = sirf(1:8:end)';
            if ~strcmp(datesall(Xaxis.freq(1),end),'1')
                disp('NOTE: Xaxis labels do not mark first quarter')
            end
            Xaxis.ticklabels = datesall(Xaxis.freq,1:4);
            Xaxis.ticklabelsize = 14;
            Yaxis.axislabelsize = 14;
            Yaxis.ticklabelsize = 10;
        else
            Title.name = varnames(V_a);
            Title.size = 11;
            Xaxis.freq = sirf(1:8:end)';
            
            if ~strcmp(datesall(Xaxis.freq(1),end),'1')
                disp('NOTE: Xaxis labels do not mark first quarter')
            end
            Xaxis.ticklabels = datesall(Xaxis.freq,1:4);
            Xaxis.ticklabelsize = 6;
            
            Title.size = 18;
            Yaxis.axislabelsize = 16;
            Yaxis.ticklabelsize = 12;
        end
        
        lgnd = [];
        
        %% Shock Decomposition
    case {'Shock Decomposition','Shock Decomposition-noDet'}
        plotSeparate = 0;
        Xaxis.limits = [min(sirf_shockdec),max(sirf_shockdec)];
        Xaxis.freq = sirf_shockdec(1:4:end);
        Xaxis.ticklabels = datesall(Xaxis.freq,1:4);
        
        if V_a==3
            Yaxis.axislabel = 'Percent';
        else
            Yaxis.axislabel = 'Percent Q-to-Q Annualized';
        end
        
        if any(V_a == Vseq_alt)
            Yaxis.axislabelsize = 12;
            Yaxis.ticklabelsize = 12;
        else
            Yaxis.axislabelsize = 11;
            Yaxis.ticklabelsize = 11;
        end
        
        if strcmp(plotType,'Shock Decomposition-noDet')
            if V_1 == length(Vseq)-length(Vseq_alt)
                lgnd.size = 11;
            elseif V_1 == length(Vseq)
                lgnd.size = 11;
            else lgnd = [];
            end
        elseif strcmp(plotType,'Shock Decomposition')
            if V_1 == length(Vseq)-length(Vseq_alt)
                lgnd.size = 14;
            elseif V_1 == length(Vseq)
                lgnd.size = 14;
            else 
                lgnd = [];
            end
        end
        
        Title.name = [varnames(V_a), '(deviations from mean)'];
        Title.size = 11;
        
        
        %% Counterfactuals
    case {'Counterfactual by Shock', 'Counterfactual by Variable'}
        plotSeparate = 0;
        % Counterforecast=1 plots 60QAhead, Counterforecast=1 or 2 plots up to stime only
        if Counterforecast==1,tick_counter=3; else tick_counter=1; end
        
        sirf_counter = (Startdate:Enddate);
        
        if strcmp(plotType, 'Counterfactual by Shock')
            if pass
                S_a = varargin{1};
                Yaxis.axislabel = varnames_YL(S_a);
                Title.name = varnames(S_a);
            end
        else
            if pass
                S_a = varargin{1};
                Title.name = shocksnames(S_a);
            end
            Yaxis.axislabel = varnames_YL(V_a);
        end
        
        Title.size = 10;
        Xaxis.limits = [min(sirf_counter),max(sirf_counter)];
        Xaxis.freq = sirf_counter(1:4*tick_counter:end);
        Yaxis.axislabelsize = 8;
        Yaxis.ticklabelsize = 6;
        % 			if plotSeparate_counter
        if plotSeparate
            Xaxis.ticklabels = datesall(Xaxis.freq,:);
        else
            Xaxis.ticklabels = datesall(Xaxis.freq,3:end);
        end
        lgnd=[];
        
        %% Shock innovations
    case {'Shock'}
        plotSeparate = 1;
        sirf_shock = (Startdate:Idate + logical(peachcount));
        
        Title.name = [shocksnames(V_a)];
        Title.size = 10;
        
        Xaxis.limits = [min(sirf_shock),max(sirf_shock)];
        
        Xaxis.freq = sirf_shock(1:4:end);
        Yaxis.axislabel = 'Standard Deviations';
        Yaxis.axislabelsize = 8;
        Yaxis.ticklabelsize = 6;
        
        % 			if plotSeparate_shocks
        if plotSeparate
            Xaxis.ticklabels = datesall(Xaxis.freq,1:4);
        else
            Xaxis.ticklabels = datesall(Xaxis.freq,3:end);
        end
        
        lgnd=[];
    case {'ShockAnt'}
        plotSeparate = 1;
        sirf_shock = (Startdate:Idate + logical(peachcount));
        
        Title.name = [shocksnames(V_a)];
        Title.size = 10;
        
        Xaxis.limits = [min(sirf_shock),max(sirf_shock)];
        
        Xaxis.freq = sirf_shock(1:4:end);
        Yaxis.axislabel = 'Percent';
        Yaxis.axislabelsize = 8;
        Yaxis.ticklabelsize = 6;
        
        % 			if plotSeparate_shocks
        if plotSeparate
            Xaxis.ticklabels = datesall(Xaxis.freq,1:4);
        else
            Xaxis.ticklabels = datesall(Xaxis.freq,3:end);
        end
        
        lgnd=[];
        
    case {'Shock Squared','Htil','Eta Squared','Sigma'}
        plotSeparate = 1;
        sirf_shock_sq = 1:nobs;
        
        Title.name = [shocksnames(V_a)];
        Title.size = 10;
        
        Xaxis.limits = [min(sirf_shock_sq),max(sirf_shock_sq)];
        
        Xaxis.freq = sirf_shock_sq(1:4:end);
        Yaxis.axislabel = 'Standard Deviations';
        Yaxis.axislabelsize = 8;
        Yaxis.ticklabelsize = 6;
        
        % 			if plotSeparate_shocks
        if plotSeparate
            Xaxis.ticklabels = datesall(Xaxis.freq,:);
        else
            Xaxis.ticklabels = datesall(Xaxis.freq,3:end);
        end
        
        lgnd=[];
        
end

line.width = 2.2;
switch plotType
    case {'Forecast' 'Forward Guidance'}
        plotSeparate = 1;
    case 'Counterfactual by Shock'
        plotSeparate = 0;
    case 'Counterfactual by Variable'
        plotSeparate = 0;
    case {'Shock','ShockAnt'}
        plotSeparate = 1;
    case {'Shock Decomposition','Shock Decomposition-noDet'}
        plotSeparate = 1;
end

