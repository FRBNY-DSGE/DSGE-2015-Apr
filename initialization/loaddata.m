%% loaddata
%
% Description: reads time series data from ASCII 
% Output: 
%        1) YYall ---> matrix of observables
%        2) XXall ---> 
%        3) ti --> 
%        4) nobs --> number of periods in data imported
%        4) dlpopall --> log differences of the population obtained from Haver Analaytics
%        5) dlMA_pop --> log differences of the population forecasted by Macro Advisers 

function [YYall,XXall,ti,nobs,dlpopall,dlMA_pop,MA_pop,population] = loaddata(nvar,nlags,nant,antlags,psize,zerobound,peachflag,mspec)

if mspec==990
    start_date = 1959.00; 
    load('rawData_990');

    popreal=data(1:end,11);
    m = length(MA_pop);
    [xxx,popfor]=Hpfilter([popreal;MA_pop(2:end)],1600);
    population = popfor(1:end-m+1);
    MA_pop = popfor(end-m+1:end);

    dlpopreal = log(popreal(2:end)) - log(popreal(1:end-1));
    dlpopall = log(population(2:end)) - log(population(1:end-1));
    dlMA_pop = log(MA_pop(2:end)) - log(MA_pop(1:end-1));
    
    nobs=size(data,1);

    ti     = (start_date:.25:(start_date+.25*(nobs-1)))';

    %% Series
    series  = nan(nobs-1,nvar);

    %Trim series for growth rate transformations
    nobs = nobs - 1;
    ti = ti(2:end);

    %% Output growth (log approximation quarterly)
    %% Load levels
    Output = 100*(log(data(2:end,1))-log(data(1:end-1,1)));
    series(:,1) = Output + 100*(dlpopreal - dlpopall); % Add back in the adjustment from HP Filtering

    %% Employment/Hours per capita
        % log hours per capita
        series(:,2) = (log(3*(data(2:end,14).*data(2:end,13)/100))-log(population(2:end)))*100;
    %% Real Wage Growth
        % per hour
        series(:,3) = 100*((log(data(2:end,15))-log(data(1:end-1,15)))-(log(data(2:end,2))-log(data(1:end-1,2))));
        nn=3;
    
    %% Inflation

    % Core PCE
    series(:,nn+1) = 100*(log(data(2:end,12))-log(data(1:end-1,12)));
        % GDP Deflator
        series(:,nn+1) = 100*(log(data(2:end,2))-log(data(1:end-1,2)));

    % Add Core PCE for models with factor structure on inflation.
    if any(mspec == [988 990 991 992])
      nn=nn+1;
      series(:,nn+1) = 100*(log(data(2:end,12))-log(data(1:end-1,12)));
    end

    %% nominal short-term interest rate (3 months) - % annualized
    %series(:,6) = robs(68:end);
    series(:,nn+2) = data(2:end,8);

    %% Consumption growth (log approximation quarterly annualized)
    %series(:,7) =  dc(68:end);
    Consumption =  100*(log(data(2:end,3))-log(data(1:end-1,3)));
    series(:,nn+3) = Consumption + 100*(dlpopreal - dlpopall);


    %% Investment growth (log approximation quarterly annualized)
    %series(:,8) =  dinve(68:end);
    Investment = 100*(log(data(2:end,4))-log(data(1:end-1,4)));
    series(:,nn+4) = Investment + 100*(dlpopreal - dlpopall);

    %% Unemployment
      nn = nn+4;


    %% spread: BAA-10yr TBill
        series(:,nn+1) = data(2:end,10); %%add for model with Financial Frictions

    %% Long Term Inflation Expectations
            series(:,nn+2) = longinf(2:end)/4;

    % Add long rate data
      data_add = [nan(length(series)-length(lrate),1); lrate]/4;
      series(:,nn+3) = data_add; % Last column is nans, just waiting for long rate data, so throw it out

    % Add output level
    if any(mspec == [977 980 981 982])
        % output log levels
        series(:,nn+4) = 100*log(data(2:end,1)) + 100*(log(popreal(2:end))-log(population(2:end)));
    end

    % Add Fernald TFP series
        Alpha=fernaldTFP(1:nobs,1);
        tfp_unadj=fernaldTFP(1:nobs,2);
        tfp_adj=fernaldTFP(1:nobs,3);
                tfp = tfp_unadj;
                  tfp=(tfp-nanmean(tfp)) ./ (4*(1-Alpha));
        series(:,nn+4) = tfp;


    %% Adding OIS expectations
        %getting expectations
        [ExpFFR,peachdata_FFR] = ExpFFR_OIS(nant,antlags,psize,zerobound,peachflag);

        % Account for fact that our dataset might stop earlier than OIS, which
        % means we should truncate OIS
        ExpFFR = ExpFFR(1:antlags+1,:);

        E_num = size(ExpFFR,1);
        fill = nan(length(series)-E_num,nant);
        data_add = [fill;ExpFFR]/4;
        series = [series,data_add];
        nvar = nvar + nant;

    YYall     = series(1+nlags:nobs,:);%%nlags switched from T0, DF
    XXall     = ones(nobs-nlags, 1+nvar*nlags);
    ti        = ti(1+nlags:nobs,:);
    dlpopall  = dlpopall(1+nlags:nobs); % NEW

    for i = 1:1:nlags;
        XXall(:,1+(i-1)*nvar+1:1+i*nvar) = series(nlags-i+1:nobs-i,:);
    end

    if nlags > 0
      XXall = [XXall;[1,YYall(end,:),XXall(end,1+1:1+nvar*nlags-nvar)]];
    end
    nobs = size(YYall,1);


    
else
    error('Appropriate data transformations not specified for given model specification. Please see loaddata.m');
    
end

