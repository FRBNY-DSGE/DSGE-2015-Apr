% Record numbers plotted
function [record] = recordFcn(plotType,record,V_1,V_a,peachcount,peachfile,...
                              tiall,Startdate,stime,qahead,...
                              varnames,names_shocks,fid0,finalRecord,varargin)
                          
    
    if V_1==1, j=1; else j=0; end % Applies to everything except Counterfactual by Shock and Counterfactual by Var
    
    switch plotType 
        
        case 'Q4Q4 Table'
            percent = varargin{1}; 
            table = varargin{2}; 
            Means_q4q4 = varargin{3}; 
            Bands_q4q4 = varargin{4}; 
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional Q4Q4 Table';
                case 0, record(end+j).type = 'Unconditional Q4Q4 Table';
                case 2, record(end+j).type = 'Semi-Conditional Q4Q4 Table';
            end
            record(end).labels = {'dates',[num2str(100*percent),' %LB'],'mean',[num2str(100*percent),' %UB']};
            record(end).data(:,:,V_1) = [table' squeeze(Bands_q4q4(1,V_1,:,peachcount+1)), Means_q4q4(V_1,:,peachcount+1)',squeeze(Bands_q4q4(2,V_1,:,peachcount+1))];
            
        case 'Forecast'
            bandsp = varargin{1}; 
            mlp = varargin{2}; 
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional Forecast';
                case 0, record(end+j).type = 'Unconditional Forecast';
                case 2, record(end+j).type = 'Semi-Conditional Forecast';
            end
            record(end).labels = {'dates','90%LB','80%LB','70%LB','60%LB','50%LB','50%UB','60%UB','70%UB','80%UB','90%UB','mean'};
            record(end).data(:,:,V_a) = [tiall(Startdate:(stime+qahead),:),bandsp(Startdate:end,:),mlp(Startdate:end)];
            
        case 'Exp_Forecast'
            bandsp = varargin{1}; 
            mlp = varargin{2}; 
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional Forecast';
                case 0, record(end+j).type = 'Unconditional Forecast';
                case 2, record(end+j).type = 'Semi-Conditional Forecast';
            end
            record(end).labels = {'dates','90%LB','80%LB','70%LB','60%LB','50%LB','50%UB','60%UB','70%UB','80%UB','90%UB','mean'};
            record(end).data(:,:,V_a) = [tiall(Startdate:(stime+qahead),:),bandsp(Startdate:end,:),mlp(Startdate:end)];
            bandsp = varargin{1}; 
            mlp = varargin{2}; 
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional Forecast';
                case 0, record(end+j).type = 'Unconditional Forecast';
                case 2, record(end+j).type = 'Semi-Conditional Forecast';
            end
            record(end).labels = {'dates','90%LB','80%LB','70%LB','60%LB','50%LB','50%UB','60%UB','70%UB','80%UB','90%UB','mean'};
            record(end).data(:,:,V_a) = [tiall(Startdate:(stime+qahead),:),bandsp(Startdate:end,:),mlp(Startdate:end)];
            
        case 'Shock Decomposition'
            list = varargin{1};
            mlp = varargin{2};
            mlp_bar = varargin{3};
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional Shock Decomposition';
                case 0, record(end+j).type = 'Unconditional Shock Decomposition';
                case 2, record(end+j).type = 'Semi-Conditional Shock Decomposition';
            end
            record(end).labels= {'dates','deviations from mean',list{:}};
            %record(end).data(:,:,V_a) = [tiall(Startdate:(stime+qahead),:),mlp(Startdate:end),mlp_bar(Startdate:end,:)];
            
        case '4Q Forecast'
            bandsp = varargin{1}; 
            mlp = varargin{2}; 
            
            switch peachcount
                case 1, record(end+j).type = 'Conditional 4Q Forecast';
                case 0, record(end+j).type = 'Unconditional 4Q Forecast';
                case 2, record(end+j).type = 'Semi-Conditional 4Q Forecast';
            end
            record(end).labels= {'dates','90%LB','80%LB','70%LB','60%LB','50%LB','50%UB','60%UB','70%UB','80%UB','90%UB','mean'};
            record(end).data(:,:,V_a) = [tiall(Startdate:(stime+qahead),:),bandsp(Startdate:end,:),mlp(Startdate:end)];
            
        case 'Counterfactual by Variable'
            S_1 = varargin{1};
            S_a = varargin{2};
            Shseq = varargin{3};
            mlp_counter_byvar = varargin{4};
            
            if V_1==1 && S_1==1, j=1; else j=0; end 
            
            switch peachcount 
                case 1, record(end+j).type = 'Conditional Counterfactual by Variable Forecast';
                case 0, record(end+j).type = 'Unconditional Counterfactual by Variable Forecast';
                case 2, record(end+j).type = 'Semi-Conditional Counterfactual by Variable Forecast';
            end
            
            shock_label = names_shocks(Shseq);
            record(end).labels= {'dates',shock_label{:}};
            if S_1==1
                record(end).data(:,1,V_a) = tiall(Startdate:(stime+qahead),:);
            end
            record(end).data(:,S_1+1,V_a) = [mlp_counter_byvar(Startdate:end)];
                        
        case 'Counterfactual by Shock'
            S_1 = varargin{1};
            S_a = varargin{2};
            Shseq = varargin{3};
            mlp_counter = varargin{4};
            
            if S_1==1 && V_1==1, j=1; else j=0; end
            
            switch peachcount 
                case 1, record(end+j).type = 'Conditional Counterfactual by Shock Forecast';
                case 0, record(end+j).type = 'Unconditional Counterfactual by Shock Forecast';
                case 2, record(end+j).type = 'Semi-Conditional Counterfactual by Shock Forecast';
            end
            
            var_label = varnames(Shseq);
            record(end).labels= {'dates',var_label{:}};
            if S_1==1
                record(end).data(:,1,V_a) = tiall(Startdate:(stime+qahead),:);
            end
            record(end).data(:,S_1+1,V_a) = [mlp_counter(Startdate:end)];
            
        case 'Shock'
            Vseq = varargin{1};
            mlp = varargin{2};
            pnobss = varargin{3};
            
            switch peachcount 
                case 1, record(end+j).type = 'Conditional Shock Innovation';
                case 0, record(end+j).type = 'Unconditional Shock Innovation';
                case 2, record(end+j).type = 'Semi-Conditional Shock Innovation';
            end
            
            shock_label = names_shocks(Vseq);
            record(end).labels= {'dates',shock_label{:}};
            if V_1==1
                record(end).data(:,1) = tiall(Startdate:(stime+(logical(peachcount)*pnobss)),:);
            end
            record(end).data(:,V_1+1) = mlp(Startdate:(stime+(logical(peachcount)*pnobss)),:);
            
    end

if finalRecord
    % Record Dick Peach's forecast
    load(peachfile)
    fprintf(fid0, 'Conditional Data \n');
    fprintf(fid0, '%s \n', peachfile);
    eval(['fprintf(fid0, ''' repmat(' %s',1,length(data)) ' \n '', varnames{:});'])
    eval(['fprintf(fid0, ''' repmat(' %4.5f',1,length(data)) ' \n '', data);'])
    fprintf(fid0, ' \n');
    
    % Record plotted data
    for k = 1:length(record)
        fprintf(fid0,'%s \n\n', record(k).type);
        for v = 1:size(record(k).data,3)
            if ~ismember(record(k).type,{'Semi-Conditional Shock Innovation';'Conditional Shock Innovation';'Unconditional Shock Innovation'})
                if ismember(record(k).type,{'Semi-Conditional Counterfactual by Shock Forecast';'Conditional Counterfactual by Shock Forecast';'Unconditional Counterfactual by Shock Forecast'})
                    fprintf(fid0,'%s \n', char(names_shocks(v)));
                else
                    fprintf(fid0,'%s \n', char(varnames(v)));
                end
            end
            fprintf(fid0,'%s ', record(k).labels{:});
            fprintf(fid0,' \n');
            for a = 1:size(record(k).data,1)
                eval(['fprintf(fid0, ''' repmat(' %4.5f',1,size(record(k).data,2)) ' \n '' ,[record(k).data(a,:,v)]);']);
            end
            fprintf(fid0,' \n \n');
        end
    end
end


end
