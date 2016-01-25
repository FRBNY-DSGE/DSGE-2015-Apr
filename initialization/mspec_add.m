
function [nvar,varnames,graph_title,cum_for,popadj,varnames_YL,varnames_irfs,varnames_YL_4Q,varnames_YL_irfs,...
    names_shocks,names_shocks_title,nl_shocks_title,shocksnames,cum_irf,vardec_varnames,shockcats,list,shockdec_color] = mspec_add(mspec,dataset,zerobound,nant,fourq_flag)

eval(['states',num2str(mspec)]);

if any(mspec==990)
    nvar = 12;
    
    % add TFP series
    if any(mspec==990)
        tfp = 'Unadjusted';
    else
        error('Must specify which TFP data series to use');
    end
    
    varnames= {'Output Growth'; 'Aggregate Hours Growth'; 'Real Wage Growth'; 'GDP Deflator';'Core PCE Inflation';...
        'Interest Rate'; 'Consumption Growth'; 'Investment Growth'; 'Spread'; 'Long Inf';...
        'Long Rate'; ['Total Factor Productivity, Util.' tfp]};
    
    varnames_irfs = { 'Output Growth'; 'Aggregate Hours Growth'; 'Real Wage Growth';...
        'GDP Deflator'; 'Core PCE Inflation'; 'Interest Rate'; 'Consumption Growth';...
        'Investment Growth'; 'Spread'; 'Long Inf'; 'Long Rate'; 'TFP'};
    
    % Save titles for graphs
    graph_title = { 'Output'; 'Labor_Supply'; 'Real_Wage'; 'GDP_Deflator'; 'Core_PCE';...
        'Interest_Rate'; 'Consumption'; 'Investment'; 'Spread'; 'Long_Inf';...
        'Long_Rate'; ['Total_Factor_Productivity_Util' tfp]};
    
    % Adjustments to be made before plotting, like going to quarterly annualized
    cum_for = [1 2 1 1 1 0 1 1 1 1 0 0];
    
    % Whether to population adjust the variables
    popadj  = [1 1 0 0 0 0 1 1 0 0 0 0];
    
    % Y-Axis Labels
    varnames_YL = {'Percent Q-to-Q Annualized'; 'Percent Q-to-Q Annualized';...
        'Percent Q-to-Q Annualized'; 'Percent Q-to-Q Annualized'; ...
        'Percent Q-to-Q Annualized';...
        'Percent Annualized'; 'Percent Q-to-Q Annualized'; 'Percent Q-to-Q Annualized';...
        'Level'; 'Level'; 'Percent Annualized'; 'Percent Annualized'};
    varnames_YL_4Q = {'Percent 4Q Growth';'Percent 4Q Growth'; 'Percent 4Q Growth';...
        'Percent 4Q Growth';'Percent 4Q Growth'; 'Percent Annualized'; 'Percent Annualized'};
    varnames_YL_irfs = {'Percent Annualized'; 'Percent Annualized'; 'Percent Annualized'; ...
        'Percent Annualized'; 'Percent Annualized'; 'Percent Annualized'; 'Percent Annualized'; ...
        'Percent Annualized'; 'Percent Annualized';...
        'Percent Annualized'; 'Percent Annualized'; 'Percent Annualized'};
    
    % Plot titles for plots of shock histories
    names_shocks = {'g'; 'b'; '\mu'; 'z'; '\lambda_f'; '\lambda_w'; 'r_m'; '\sigma_w'; '\mu_e';...
        '\gamma'; '\pi_*'; 'e_{LR}'; 'z^p'; 'e_{TFP}'; 'e_{gdpdef}'; 'e_{pce}'};...
        names_shocks_title = {'govt'; 'asset'; 'inv_tech'; 'neutral_tech'; 'price_mkp'; 'wage_mkp';...
        'Money'; 'spread'; 'mue'; 'gamma'; 'pistar'; 'me_LR'; 'zp'; 'me_tfp'; ...
        'me_gdpdef'; 'me_pce'};
    shocksnames = names_shocks';
    
    % Save names for plots of shock histories
    nl_shocks_title = {'Technology'; 'phi'; 'Financial'; 'b'; 'g'; 'Mark-Up'; 'Monetary Policy'};
    
    
    % Not relevant
    cum_irf = [];
    vardec_varnames = {'Output Growth'; 'Aggregate Hours Growth'; 'Real Wage';...
        'GDP Deflator'; 'Core PCE Inflation'; 'Interest Rate'; 'Consumption Growth';...
        'Investment Growth';};
    
    rm_total_sh=rm_sh;
    if nant>0
        for i = 1:nant
            eval(strcat('rm_total_sh = [rm_total_sh rm_shl',num2str(i),'];'));
        end
    end
    
    % Which shocks to plot in the shock decompositions:
    % - Give the indices of the shocks, where the labels (used below) aliasing
    %   the index numbers come from running statesMSPEC above
    % - Can group shocks by making an entry be an array like [mu_sh z_sh]
    shockcats = {g_sh;b_sh;[ gamm_sh mue_sh sigw_sh]; z_sh; laf_sh; law_sh; rm_total_sh; ...
        pist_sh; mu_sh; [lr_sh tfp_sh gdpdef_sh pce_sh]; zp_sh; 0;};
    
    % How to label each element of shockcats in the legend
    list = {'g'; 'b'; 'FF'; 'z'; 'p-mkp'; 'w-mkp'; 'pol.'; 'pi-LR'; 'mu'; 'me'; 'zp'; 'dt';};
    
    % Colors of the bars in the shockdecs
    shockdec_color = {'firebrick'; [0.3 0.3 1.0]; 'indigo'; 'darkorange'; 'yellowgreen'; 'teal';...
        'gold'; 'pink'; 'cyan'; [0 0.8 0]; [0 0.3 0]; 'gray';}';
    
    
    %% Adding information for expectations
    
    if zerobound
        gen_temp = @(pattern) arrayfun(@(i_) sprintf(pattern, i_), 1:nant, 'UniformOutput', false);
        temp1 = gen_temp('exp_%i');
        temp2 = gen_temp('');
        temp3 = gen_temp('Ant %i');
        temp4 = gen_temp('%iqA ant.sh.');
        temp5 = gen_temp('Ant%i');
        
        varnames(end+1:end+nant)           = temp1(1:nant);
        varnames_YL(end+1:end+nant)        = temp2(1:nant);
        varnames_YL_4Q(end+1:end+nant)     = temp2(1:nant);
        varnames_irfs(end+1:end+nant)      = temp2(1:nant);
        varnames_YL_irfs(end+1:end+nant)   = temp2(1:nant);
        cum_for(end+1:end+nant)            = zeros(1,nant);
        popadj(end+1:end+nant)             = zeros(1,nant);
        names_shocks(end+1:end+nant)       = temp3(1:nant);
        nl_shocks_title(end+1:end+nant)    = temp4(1:nant);
        names_shocks_title(end+1:end+nant) = temp5(1:nant);
        
        shocksnames = names_shocks';
    end
end
