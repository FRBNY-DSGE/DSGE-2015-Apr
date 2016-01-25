function yp = getyp(ypath,qahead,I,cum_for,popadj,dlpop_fcst,YY,Idate,q_adj,varargin)

if nargin == 10
    pop_fcst = varargin{1};
end

yp = ypath(:,qahead*(I-1)+1:qahead*I);
if cum_for(I) == 1 && popadj(I) == 0
    % Transform from log growth rates to % growth rates (annualized)
    
    % This should only be used for output, consumption, investment
    % and GDP deflator (inflation).  If you have any doubt, please
    % check the loaddata file.
    %yp = 100*4*(exp(yp/100)-1);
    
    yp = 100*((exp(yp/q_adj)).^4-1);

elseif cum_for(I) == 1 && popadj(I) == 1
    %transform from log growth rates to total (not per-capita) % growth
    %rates (annualized).
    %  yp = 100*(exp(yp/100 + .01106) - 1); % population estimate
    %  yp = 100*(exp(yp/q_adj + dlMA_pop(end)).^4 - 1); %
    %  population estimate
    
    yp = 100*(exp(yp/q_adj + dlpop_fcst).^4 - 1); % population estimate
        
    %yp=1/400*yp;
    %yp=exp(yp);
    %yp=100*(yp-1);
    %yp=yp.^4+1.106;

elseif cum_for(I) ==2 && popadj(I) == 0
    % Transform from log level to 4-quarter annualized percent change
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
    yp = ((exp(yp(:,1:end)/100 - [repmat(YY(Idate,I), [size(yp,1) 1]) yp(:,1:end-1)]/100).^4)-1)*100;

elseif cum_for(I) == 2 && popadj(I) == 1

    yp = ((exp(yp(:,1:end)/100 - [repmat(YY(Idate,I), [size(yp,1) 1]) yp(:,1:end-1)]/100 + dlpop_fcst).^4)-1)*100;

elseif cum_for(I)==3
    % Transform from log levels to level
    % The lines below specifically applies to labor share
    % Refer to loaddata for data transformation.
    % Be careful of using this for any other observables.
    yp = exp(yp/100);
elseif cum_for(I) == 4
    % This transforms log perecentage levels back to 
    % percentage levels.  This is specifically for unemployment
    % which enters the model as log-unemployment
    yp=(exp(yp/100)-1)*100;
elseif cum_for(I) == 5
    yp = yp*4;
    
elseif cum_for(I) ==6 && popadj(I) ==1
    yp = 100*(yp/q_adj + dlpop_fcst);
elseif cum_for(I) ==6 && popadj(I) ==0
    yp = 100*(yp/q_adj);

elseif cum_for(I) == 7
    % Output is in Log-levels
    yp = yp;
    
elseif cum_for(I) ==0 && popadj(I) ==1
   yp = yp + 100*log(pop_fcst);
else
    yp = yp*400/q_adj;
end;
