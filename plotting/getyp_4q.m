function yp = getyp_4q(ypath,qahead,I,cum_for,popadj,dlpop_fcst,YY,Idate,dlpop_old,q_adj)

yp = ypath(:,qahead*(I-1)+1:qahead*I);
if cum_for(I) == 1 && popadj(I) == 0
    % Transform from log growth rates to 4-Quarter percent growth rates
    % This should only be used for output, consumption, investment
    % and GDP deflator (inflation).  If you have any doubt, please
    % check the loaddata file.

    % More data is needed for 4Q growth rates.
    yp_all = [repmat(YY(end-2:end,I)',size(yp,1),1) yp];
    yp_4q = yp_all(:,end-qahead+1:end) + yp_all(:,end-qahead:end-1) + yp_all(:,end-qahead-1:end-2) + yp_all(:,end-qahead-2:end-3);
    yp = 100*(exp(yp_4q/q_adj) - 1);

elseif cum_for(I) == 1 && popadj(I) == 1
    %transform from log growth rates to total (not per-capita) % growth
    %rates (annualized)

    yp_all = [repmat(YY(end-2:end,I)',size(yp,1),1) yp];
    yp_4q = yp_all(:,end-qahead+1:end) + yp_all(:,end-qahead:end-1) + yp_all(:,end-qahead-1:end-2) + yp_all(:,end-qahead-2:end-3);
    
    dlpop_all = [dlpop_old(:,end-2:end) dlpop_fcst];
    dlpop_4q = dlpop_all(:,end-qahead+1:end) + dlpop_all(:,end-qahead:end-1) + dlpop_all(:,end-qahead-1:end-2) + dlpop_all(:,end-qahead-2:end-3);
    yp = 100*(exp(yp_4q/q_adj + dlpop_4q) - 1);

elseif cum_for(I) ==2 && popadj(I) == 0
    % Transform from log level to 4-Quarter percent change
    % The lines below specifically applies to labor supply (hours worked
    % per capita). Refer to loaddata for the data transformation.
    % Be careful of using this for any other observables.
    % Some notes on what is happening in the line directly below:
    % YY is necessary to get the last data point so that a percent
    % change can be computed for the first period. Repmat is used
    % to put the data point in each row of the simulations.  The
    % log levels are subtracted to get the log percent changes and
    % then the exponential is used to remove the log from the
    % levels.
    yp_old = [repmat(YY(end-3:end,I)',size(yp,1),1) yp(:,1:end-4)];
    yp_4q = yp - yp_old;
    yp = 100*(exp(yp_4q/100)-1);
%     yp = ((exp(yp(:,1:end)/100 - [repmat(YY(Idate,I), [size(yp,1) 1]) yp(:,1:end-1)]/100).^4)-1)*100;

elseif cum_for(I) == 2 && popadj(I) == 1

%     yp = ((exp(yp(:,1:end)/100 - [repmat(YY(Idate,I), [size(yp,1) 1]) yp(:,1:end-1)]/100 + dlpop_fcst).^4)-1)*100;

    yp_old = [repmat(YY(end-3:end,I)',size(yp,1),1) yp(:,1:end-4)];
    yp_4q = yp - yp_old;
    
    dlpop_all = [dlpop_old(:,end-2:end) dlpop_fcst];
    dlpop_4q = dlpop_all(:,end-qahead+1:end) + dlpop_all(:,end-qahead:end-1) + dlpop_all(:,end-qahead-1:end-2) + dlpop_all(:,end-qahead-2:end-3);
    
    yp = 100*(exp(yp_4q/100 + dlpop_4q) - 1);

elseif cum_for(I)==3
    % Transform from log levels to level
    % The lines below specifically applies to labor share
    % Refer to loaddata for data transformation.
    % Be careful of using this for any other observables.
    yp = exp(yp/100);
    
elseif cum_for(I) == 4
    
    yp_all = [repmat(YY(end-2:end,I)',size(yp,1),1) yp];
    yp = (yp_all(:,end-qahead+1:end) + yp_all(:,end-qahead:end-1) + yp_all(:,end-qahead-1:end-2) + yp_all(:,end-qahead-2:end-3))/4;
elseif cum_for(I) == 5
    yp = yp*4;
else
    
    yp = yp*400/q_adj;
    
end;