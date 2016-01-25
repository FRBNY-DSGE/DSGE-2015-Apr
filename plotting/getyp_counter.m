function yp_counter = getyp_counter(ypath_counter,counter_ahead,I,cum_for,popadj,mnobss,dlpop_old,dlpop_ctr,yend,q_adj)

% This term should be 400 assuming that growth rates are log quarterly
% growth rates multiplied by 100
%if q_adj ~= 400, warning('q_adj is not 400'); end

yp_counter = ypath_counter(:,counter_ahead*(I-1)+1:counter_ahead*I);
if cum_for(I) == 1 && popadj(I) == 0
    % Transform from log growth rates to % growth rates (annualized)
    %yp = 100*4*(exp(yp/100)-1);
    yp_counter = 100*((exp(yp_counter/q_adj)).^4-1);
elseif cum_for(I) == 1 && popadj(I) == 1
    % Transform from log growth rates to total (not per-capita) % growth
    % rates (annualized)
%                     yp_counter(:,mnobss+1:mnobss+size(peachdata,1)) = 100*(exp(yp_counter(:,mnobss+1:mnobss+size(peachdata,1))/q_adj + dlMA_pop).^4 - 1); % population estimate
%                     yp_counter(:,1:mnobss) = 100*(exp(yp_counter(:,1:mnobss)/400 + repmat(dlpop(end-mnobss+1:end)',5000,1)).^4 - 1);

    yp_counter(:,mnobss+1:end) = 100*(exp(yp_counter(:,mnobss+1:end)/q_adj + dlpop_ctr).^4 - 1); % population estimate
    yp_counter(:,1:mnobss)= 100*(exp(yp_counter(:,1:mnobss)/q_adj + dlpop_old).^4 - 1);
elseif cum_for(I)==2 && popadj(I) == 0
    % Transform from log level to 4-quarter annualized percent change
    % The lines below specifically applies to labor supply (hours worked
    % per capita). Refer to loaddata for data transformation.
    % Be careful of using this for any other observables.
    % Some notes on what is happening in the line directly below:
    % YY is necessary to get the last data point so that a percent
    % change can be computed for the first period. Repmat is used
    % to put the data point in each row of the simulations.  The
    % log levels are subtracted to get the log percent changes and
    % then the exponential is used to remove the log from the
    % levels.
%         yp_counter = ((exp(yp_counter(:,1:end)/100 - [repmat(YY(Startdate,I), [size(yp_counter,1) 1]) yp_counter(:,1:end-1)]/100).^4)-1)*100;
    yp_counter = ((exp(yp_counter(:,1:end)/100 - [yend(:,I) yp_counter(:,1:end-1)]/100).^4)-1)*100;
        %yp_counter = ((exp(yp_counter(:,1:end)/100 - [repmat(YY(Idate,I), [size(yp_counter,1) 1]) yp_counter(:,1:end-1)]/100).^4)-1)*100;
elseif cum_for(I) == 2 && popadj(I) == 1
    yp_counter(:,mnobss+1:end) = ((exp(yp_counter(:,mnobss+1:end)/100 - yp_counter(:,mnobss:end-1)/100 + dlpop_ctr).^4)-1)*100;
    yp_counter(:,1:mnobss) = ((exp(yp_counter(:,1:mnobss)/100 - [yend(:,I) yp_counter(:,1:mnobss-1)]/100 + dlpop_old).^4)-1)*100;
elseif cum_for(I)==3
    % Transform from log levels to level
    % The lines below specifically applies to labor share
    % Refer to loaddata for data transformation.
    % Be careful if using this for any other observables.
    yp_counter = exp(yp_counter/100);
elseif cum_for(I) == 4
    % This transforms log perecentage levels back to 
    % percentage levels.  This is specifically for unemployment
    % which enters the model as log-unemployment
    yp_counter=(exp(yp_counter/100)-1)*100;
elseif cum_for(I) == 5
    yp_counter = yp_counter*4;
elseif cum_for(I) ==6 && popadj(I) ==1
    yp_counter(:,mnobss+1:end) = 100*(yp_counter(:,mnobss+1:end)/q_adj + dlpop_ctr);
    yp_counter(:,1:mnobss) = 100*(yp_counter(:,1:mnobss)/q_adj + dlpop_old);
    
elseif cum_for(I) ==6 && popadj(I) ==0
    yp_counter(:,mnobss+1:end) = 100*(yp_counter(:,mnobss+1:end)/q_adj);
    yp_counter(:,1:mnobss) = 100*(yp_counter(:,1:mnobss)/q_adj);
elseif cum_for(I) ==7 
    yp_counter(:,1:mnobss) = yp_counter(:,1:mnobss);
else
    yp_counter = yp_counter*400/q_adj;
end;
