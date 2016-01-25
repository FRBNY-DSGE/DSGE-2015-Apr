%                                   OVERVIEW
% forecastFcn_est_ant.m: a function called on by forecast_mode_est_ant.m. This
%                function takes draws allocated each worker to calculate
%                forecasts, counterfactuals, shock decompositions, and
%                smoothed
%                shock estimates, both unconditional and conditional on
%                Dick
%                Peach's forecast.
%
% forecastFcn_est_ant.m is broken up into 4 main steps:
%
% STEP 1: Get matrices of the state transition and
%         measurement equations. In gibb.m, these
%         matrices are computed and condensed, so here we read in and
%         reshape the input data.
% STEP 2: Run the Kalman filter.
% STEP 3: Run the Kalman smoother to get a sequence of smoothed states
%         (S_0/T,...,S_T/T). In this step, we also get a sequence of
%         smoothed
%         shocks, which we plot (these figures are referred to as "shock
%         histories").
% STEP 4: Using S_{T/T}, compute forecasts. With S_{0/T}, compute the
%         deterministic trend. Using a subset of the sequence of smoothed ...
%         states, we can compute the counterfactuals and shock decompositions.
%
% The state transition equation is: S_t = TTT*S_t-1+RRR*eps_t
% The measurement equation is: Y_t = ZZ*S_t+DD
%
%                           IMPORTANT INPUTS
% k: worker number index (used to ensure workers receive independently
%     randomized seeds)
%
% priotheta (npara x nsim): parameter draws
% TTTsim (nstate^2 x nsim), RRRsim (nstate*nshocks x nsim): state transition
%     equation matrices
% zendsim (nstate x nsim): end-of-sample smoothed state vector (from gibb.m)
%
% zerobound: flag specifying zerobound model (incorporates anticipated shocks)
% bdd_int_rate: flag specifying bounded interest rate rule
% nant,antlags: number of anticipated shocks, number of past periods in
%               which ZB was in effect.
% nshocks_all: total number of shocks (incl. anticipated shocks)
%
% peachflag: (if 0) no conditional data ("unconditional forecast")
%            (if 1) use conditional forecasts and data for all observables
%                   except labor share ("conditional forecast")
%            (if 2) use conditional data for spreads and FFR only
%                   ("semi-conditional forecast")
% peachdata, semi_peachdata (both are nant x (nvar+nant)): conditional,
%     semi-conditional data
%
% nobs: length of observable time series
% stime: index for the current quarter, relative to the observable time series
%
% Startdate: Index for counterfactual start date. Also, index for shock
%            decomposition figure start date.
% Enddate_forecastfile: total number of quarters in the observable time
%                       series and the forecast period: stime + qahead

%                               OUTPUTS
% For output variable dimensions, see section "Initialize Output Variables"
% below.
%
% fcast: forecast; there are subfields for unconditional, semi-conditional,
%        conditional, and (if zerobound) their bounded forecast counterparts
% varargout (optional outputs):
% shocks: smoothed shock histories
% ytrend: steady state, for shock decomposition plots
% dettrend: deterministic trend, for shock decomposition plots
% shockdec_all: shock contributions to forecast, for shock decomposition plots
% counter_all: counterfactual forecasts
%
%                           IMPORTANT SUBFUNCTIONS
% augmentFilter.m, augmentSmoother.m: augments relevant matrices for the Kalman
%     Filter and Smoother
% kalcvf2NaN.m: runs the Kalman Filter
% kalsmth_k93.m: runs the Kalman Smoother using algorithm in (Koopman 1993)


function [fcast,varargout] = forecastFcn_est_ant(k,priotheta,TTTsim,RRRsim,zendsim,...
    mspec,zerobound,bdd_int_rate,peachflag,plotImp,nimpVar,...
    nvar,qahead,qahead_extra,nshocks_all,YY,nobs,peachdata,psize,...
    Enddate_forecastfile,Startdate,stime,ShockremoveList,...
    nant,antlags,varnames,names_shocks,...
    nstate_0,nshocks_0,sflag,YY0,...
    valid_para,nlags,npara,coint,cointadd,ExpFFR,nplotstates,varargin)

% Set new seed (random number sequence) for each worker
rand('state',sum(100*clock+k))
randn('state',sum(100*clock+k))


% Number of draws distributed. Matrices holding draws have already been indexed
% by jstep in forecast_parallel.m.
nDraws = size(TTTsim,1);

%% Initialize output variables
fcast.uncond_hist = zeros(nDraws,nvar*stime);
fcast.uncond = zeros(nDraws,nvar*qahead);

if peachflag
    fcast.cond_hist = zeros(nDraws,nvar*stime);
    fcast.cond = zeros(nDraws,nvar*qahead);
end

if zerobound
    fcast.uncond_bdd = zeros(nDraws,nvar*qahead);
    if peachflag, fcast.cond_bdd = zeros(nDraws,nvar*qahead); end
end

shocks.data = zeros(nDraws,nshocks_all*size(YY',2));

% Non-standardized shocks
shocks.data_ns = zeros(nDraws,nshocks_all*size(YY',2));
if peachflag
    shocks.peach = zeros(nDraws,nshocks_all*(size(peachdata,1)));
    shocks.datawpeach = zeros(nDraws,nshocks_all*size(YY,1));
    shocks.peach_ns = zeros(nDraws,nshocks_all*(size(peachdata,1)));
    shocks.datawpeach_ns = zeros(nDraws,nshocks_all*size(YY,1));
end


ytrend = zeros(nDraws,nvar);
dettrend = zeros(nDraws,nvar*(Enddate_forecastfile - Startdate + 1),...
    peachflag+1);
counter_all = zeros(nDraws,nvar*(Enddate_forecastfile - Startdate),...
    peachflag+1,length(ShockremoveList));
shockdec_all = zeros(nDraws,nvar*(Enddate_forecastfile - Startdate + 1),...
    peachflag+1,length(ShockremoveList) - 1);


% Load the 2part classes and set up some important information
adj = adj2part(mspec);

%% Loop through parameter draws
for a = 1:nDraws
    
    nstate = nstate_0;
    nshocks = nshocks_0;
    
    % state selector
    params = priotheta(a,:)';
    [A, C_ss, cum_for, states_names] = mapStates(mspec,nstate,nplotstates,params);
    
    
    %% STEP 1: Reshape inputs to get matrices of the state transition and
    %% measurement equations
    % State transition: S_t = TTT*S_(t-1) + RRR*eps_t
    % Measurement equation: Y_t = ZZ*S_t + DD
    
    if zerobound && is2part(mspec)
        [TTT,RRR,zend] = getmodel(TTTsim,RRRsim,zendsim,a,nstate+nant+adj,nshocks+nant);
    else
        [TTT,RRR,zend] = getmodel(TTTsim,RRRsim,zendsim,a,nstate,nshocks);
    end
    
    nstate = size(TTT,1);
    
    if is2part(mspec)
      if zerobound

        % Get the starting indices
        [start_ant_state, start_ant_shock, revol_ind] = get_start_ant(mspec, nant);

        % Get noZB matrices
        [TTT_noZB, RRR_noZB, CCC_noZB]  = ...
          dsgelh_getNoZB(nant, start_ant_state, start_ant_shock, TTT, RRR, zeros(size(TTT,1),1),revol_ind{:});

        zend_noZB = zend;
        zend_noZB( start_ant_state:(start_ant_state+nant-1) ) = []; % Ditch ant states
        if ~isempty(revol_ind)
          zend_noZB(revol_ind{:}) = []; % Ditch extra equation
        end

      else
          TTT_noZB = TTT;
          RRR_noZB = RRR;
          zend_noZB = zend;
      end

      nant_rm = nant;

    else
      TTT_noZB = TTT;
      RRR_noZB = RRR;
      zend_noZB = zend;
      nant_rm = 0;
    end
    
    nstate = size(TTT_noZB,1);
    [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = getmeasur(mspec,TTT_noZB,RRR_noZB,valid_para,params,nvar,nlags,npara,coint,cointadd);

    
    %% STEP 2 (Unconditional data): Kalman Filter
    
    % If we are only doing an unconditional forecast we can just skip
    % to the end since all we need is zend, which we already loaded in.
    
    
    HH = EE+MM*QQ*MM';
    VV = QQ*MM';
    VVall = [[RRR_noZB*QQ*RRR_noZB',RRR_noZB*VV];[VV'*RRR_noZB',HH]];
    lead = 1;
    
    % Define the initial mean and variance for the state vector
    % (deals with nonstationarity on model by model basis)
    [A0,P0] = lyap_nonstationary(mspec,params,TTT_noZB,RRR_noZB,QQ);
    
    % Kalman filtering over the presample
    if ~isempty(YY0)
        [pyt0,zend,pend,pred0,vpred0] = kalcvf2NaN(YY0(:,1:nvar-nant_rm)',1,zeros(nstate,1),TTT_noZB,DD,ZZ,VVall,A0,P0);
    else
        zend = A0;
        pend = P0;
        pred0 = [];
        vpred0 = [];
    end
    
    % Kalman filtering over the main sample.
    if zerobound == 1
        
        % First, filter using normal model up to period T-antlags-1.
        [L,zend,pend,pred,vpred] = kalcvf2NaN(YY(1:end-antlags-1,1:nvar-nant_rm)',lead,zeros(nstate,1),TTT_noZB,DD,ZZ,VVall,zend,pend);
        
        
        % Now change models to incorporate anticipated shocks.
        valid_a = valid_para;
        [TTT,RRR,CCC,valid_a] = dsgesolv(mspec,params,nant);
        
        [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = getmeasur(mspec,TTT,RRR,valid_a,params,nvar,nlags,npara,coint,cointadd,nant);
        
        % Might need to check that this carries through to
        % augmentStates, where there's an if statement
        r_tl1 = getState(mspec,nant,'rm_tl1');
        r_tlx = r_tl1+1;
        
        nstate = length(TTT);
        nshocks = size(RRR,2);
        
        % state selector
        [A, C_ss] = mapStates(mspec,nstate,nplotstates,params);
        
        ind_r = find(strcmp(varnames,'Interest Rate'));
        [zprev,pprev,MM_ant,EE_ant,ZZ_e,DD_e,YY_ant] = ...
            augmentFilter(r_tl1,r_tlx,nant,antlags,ind_r,zend,pend,...
            MM(1:end-nant_rm,:),...
            EE(1:end-nant_rm,1:end-nant_rm),...
            ZZ(1:end-nant_rm,:),DD(1:end-nant_rm),...
            YY(:,1:end-nant_rm),TTT,ExpFFR);
        
        HH_ant = EE_ant+MM_ant*QQ*MM_ant';
        VV_ant = QQ*MM_ant';
        VVall_ant = [[RRR*QQ*RRR',RRR*VV_ant];[VV_ant'*RRR',HH_ant]];
        
        
        % Next, filter using ZB model for periods T-antlags:T
        % These new filtered states include anticipated shocks and
        % are consistent with the zero bound.
        [L_ant,zend,pend,pred_ant,vpred_ant] = ...
            kalcvf2NaN(YY_ant',lead,zeros(nstate,1),TTT,DD_e,ZZ_e,...
            VVall_ant,zprev,pprev);
        
    else
        % Filter over the main sample.
        [L,zend,pend,pred,vpred] = ...
            kalcvf2NaN(YY',lead,zeros(nstate,1),TTT,DD,ZZ,VVall,zend,pend);
        
    end %End of if zerobound conditional.
    
    zend_uncond = zend;
    
    
    %% STEP 3 (Unconditional data): Kalman Smoother
    
    if zerobound ==1
        [pred_all,vpred_all,YY_all,A0_,P0_] = ...
            augmentSmoother(nstate,r_tlx,r_tl1,nant,antlags,pred0,pred,...
            pred_ant,vpred0,vpred,vpred_ant,YY0,YY,YY_ant,A0,P0);
    else
        pred_all = [pred0,pred];
        
        vpred_all = zeros(nstate,nstate,size(YY0,1)+size(YY,1));
        if ~isempty(YY0)
            vpred_all(:,:,1:size(YY0,1)) = vpred0;
        end
        vpred_all(:,:,size(YY0,1)+1:end) = vpred;
        
        YY_all = [YY0;YY];
        
        A0_ = A0;
        P0_ = P0;
        
        ZZ_e = ZZ;
        DD_e = DD;
    end
    
    % This section does the smoothing. For the zero bound,
    % this can be done in one step, since there is code in
    % both the Kalman Smoother and Simulation Smoother that
    % set the QQ matrix to zero in all periods before the
    % model switch, which ensures that the smoothed shocks
    % include no anticipated shocks prior to that time.
    [sm_states_all,sm_shocks_all] = ...
        kalsmth_k93(A0_,P0_,YY_all',pred_all,vpred_all,TTT,RRR,...
        QQ,ZZ_e,DD_e,nant,antlags,0,psize);
    
    % Since we smoothed back farther than we needed to
    % (even before the presample), we only need parts of
    % sm_states_all and sm_shocks_all.
    
    z0T = zeros(nstate,peachflag+1);
    if size(YY0,1) ~= 0
        z0T(:,1) = sm_states_all(:,size(YY0,1));
    else
        z0T(:,1) = A0_;
    end
    
    % z0T is the smoothed time 0 state vector. If peachflag
    % is on it will have two columns to hold one vector
    % smoothed without peachdata, and one smoothed with
    % peachdata.
    
    sm_states = sm_states_all(:,size(YY0,1)+1:end);
    sm_shocks = sm_shocks_all(:,size(YY0,1)+1:end);
    
    % Calculate shock histories
    QQinv = zeros(nshocks);
    QQinv(QQ > 0) = QQ(QQ > 0).^(-1/2);
    sm_sdz_shocks = QQinv*sm_shocks;
    
    shocks.data_ns(a,:) = sm_shocks(:)';
    shocks.data(a,:) = sm_sdz_shocks(:)';
    
    %% STEP 2 (Conditional Data): Kalman filter
    % Since the new smoothers and filter can automatically deal
    % with NaNs in the data matrix, the code can be much
    % simplified.
    
    VVall = zeros(nstate+length(peachdata));
    VVall(1:nstate,1:nstate) = RRR*QQ*RRR';
    
    [L,zend_peach,pend_peach,pred_peach,vpred_peach] = ...
        kalcvf2NaN(peachdata',lead,zeros(nstate,1),TTT,DD_e,ZZ_e,VVall,...
        zend,pend);
    
    %% STEP 3 (Conditional Data): Kalman smoother
    
    YY_all_peach = [YY_all;peachdata];
    
    pred_all_peach = [pred_all,pred_peach];
    
    vpred_all_peach = zeros(length(TTT),length(TTT),...
        size(YY_all_peach,1));
    vpred_all_peach(:,:,1:size(YY_all,1)) = vpred_all;
    vpred_all_peach(:,:,size(YY_all,1)+1:end) = vpred_peach;
    
    
    [sm_states_all_peach,sm_shocks_all_peach] = ...
        kalsmth_k93(A0_,P0_,YY_all_peach',pred_all_peach,...
        vpred_all_peach,TTT,RRR,QQ,ZZ_e,DD_e,nant,...
        antlags,1,psize);
    
    if size(YY0,1) ~= 0
        z0T(:,2) = sm_states_all_peach(:,size(YY0,1));
    else
        z0T(:,2) = A0_;
    end
    
    sm_states_peach = sm_states_all_peach(:,size(YY0,1)+1:end);
    sm_shocks_peach = sm_shocks_all_peach(:,size(YY0,1)+1:end);
    
    peachdataimplied = ...
        getpeachdataimplied(size(peachdata,1),ZZ,DD,...
        sm_states_peach(:,end-size(peachdata,1)+1:end));
    
    % Calculate conditional shocks histories
    sm_sdz_shocks_peach = QQinv*sm_shocks_peach;
    
    shocks.datawpeach_ns(a,:) = ...
        reshape(sm_shocks_peach(:,1:size(YY,1)),1,nshocks*size(YY,1));
    shocks.peach_ns(a,:) = ...
        reshape(sm_shocks_peach(:,size(YY,1)+1:end),1,nshocks*psize);
    
    shocks.datawpeach(a,:) = ...
        reshape(sm_sdz_shocks_peach(:,1:size(YY,1)),1,nshocks*size(YY,1));
    shocks.peach(a,:) = ...
        reshape(sm_sdz_shocks_peach(:,size(YY,1)+1:end),1,nshocks*psize);
    
    %% STEP 4: Compute forecasts
    
    ind_r = find(strcmp(varnames,'Interest Rate'));
    ind_r_sh = find(strcmp(names_shocks,'r_m'));
    
    % Conditional forecast
    if peachflag
        
        if zerobound
            
            % Calculate forecast with bounded interest rate off
            [yypred,yypred_s] = getForecast_s(qahead_extra,nvar,nshocks,zeros(size(TTT,1),1),TTT,RRR,DD,ZZ,QQ,zend_peach,...
                ind_r,ind_r_sh,zerobound,0,nant,sflag,A,C_ss,nplotstates,mspec);
            
            yypred_nocoint = [peachdataimplied; ...
                yypred(:,1:nvar)];
            fcast.cond(a,:) = yypred_nocoint(:)';
            
            % Calculate forecast with bounded interest rate on
            [yypred,yypred_s] = getForecast_s(qahead_extra,nvar,nshocks,zeros(size(TTT,1),1),TTT,RRR,DD,ZZ,QQ,zend_peach,...
                ind_r,ind_r_sh,zerobound,1,nant,sflag,A,C_ss,nplotstates,mspec);
            yypred_nocoint = [peachdataimplied; ...
                yypred(:,1:nvar)];
            fcast.cond_bdd(a,:) = yypred_nocoint(:)';
            
        else
            
            [yypred] = ...
                getForecast_s(qahead_extra,nvar,nshocks,zeros(size(TTT,1),1),...
                TTT,RRR,DD,ZZ,QQ,zend_peach,ind_r,ind_r_sh,...
                zerobound,bdd_int_rate,nant,sflag,A,...
                nplotstates,mspec);
            
            yypred_nocoint = [peachdataimplied; ...
                yypred(:,1:nvar)];
            fcast.cond(a,:) = yypred_nocoint(:)';
            
            
        end %End of if zerobound conditional.
        
    end %End of if peachflag conditional.
    
    % Unconditional forecast
    if zerobound
        
        [yypred,yypred_s] = getForecast_s(qahead,nvar,nshocks,zeros(size(TTT,1),1),TTT,RRR,DD,ZZ,QQ,zend_uncond,...
            ind_r,ind_r_sh,zerobound,0,nant,sflag,A,C_ss,nplotstates,mspec);
        fcast.uncond(a,:) = yypred(:)';
        
        % Calculate bounded interest rate forecast
        [yypred,yypred_s] = getForecast_s(qahead,nvar,nshocks,zeros(size(TTT,1),1),TTT,RRR,DD,ZZ,QQ,zend_uncond,...
            ind_r,ind_r_sh,zerobound,1,nant,sflag,A,C_ss,nplotstates,mspec);
        fcast.uncond_bdd(a,:) = yypred(:)';
        
    else
        [yypred] = ...
            getForecast_s(qahead,nvar,nshocks,zeros(size(TTT,1),1),TTT,RRR,...
            DD,ZZ,QQ,zend_uncond,ind_r,ind_r_sh,zerobound,...
            bdd_int_rate,nant,sflag,A,nplotstates,mspec);
        
        fcast.uncond(a,:) = yypred(:)';
        
        
    end %End of if zerobound conditional.
    
    % Compute Counterfactuals, shock decompositions, deterministic trends, and steady states
    
    ytrend(a,:) = DD';
    
    %deterministic trends
    [dettrend(a,:,:)] =...
        getdettrend_s(Enddate_forecastfile,nvar,ZZ,TTT,DD,z0T,Startdate,peachflag,A,nplotstates);
    %counterfactuals
    [counter_all(a,:,:,:)] = ...
        getcounter_s(ShockremoveList,peachflag,psize,...
        Enddate_forecastfile,sm_states,sm_shocks,...
        sm_states_peach,sm_shocks_peach,...
        Startdate,nobs,nshocks,qahead,...
        nvar,TTT,RRR,DD,ZZ,ind_r,ind_r_sh,bdd_int_rate,A,...
        nplotstates,mspec);
    %shock decompositions
    [shockdec_all(a,:,:,:)] = ...
        getshockdec_s(ShockremoveList,peachflag,psize,...
        Enddate_forecastfile,sm_shocks,sm_shocks_peach,...
        Startdate,nobs,qahead,nvar,TTT,...
        RRR,DD,ZZ,ind_r,ind_r_sh,bdd_int_rate,A,nplotstates);
    
end

varargout{1} = shocks;
varargout{2} = ytrend;
varargout{3} = dettrend;
varargout{4} = counter_all;
varargout{5} = shockdec_all;
