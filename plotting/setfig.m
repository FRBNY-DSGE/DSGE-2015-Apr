function [figmain] = setfig(plotType)
    switch plotType
        case 'Forecast'
            figmain = 1000;
        case 'Shock Decomposition'
            figmain = 2000;
        case 'Forecast Comparison'
            figmain = 3000;
        case '4Q Forecast'
            figmain = 4000;
        case 'Counterfactual by Shock'
            figmain = 5000;
        case 'Counterfactual by Variable'
            figmain = 6000;
        case 'Shock'
            figmain = 7000;
        case 'Exp_Forecast'
            figmain = 8000;
        case 'Shock Decomposition-noDet'
            figmain = 9000;
        case 'Shock Squared'
            figmain = 10000;
        case 'Htil'
            figmain = 11000;
        case 'Eta Squared'
            figmain = 12000;
        case 'Sigma'
            figmain = 13000;
        case 'Forecast Comparison Separate'
            figmain = 15000;
        case '4Q Forecast Separate'
            figmain = 16000;
        case 'ShockAnt'
            figmain = 17000;
        case 'Forward Guidance'
            figmain = 18000;
        case 'Forward Guidance Horizontal'
            figmain = 19000;
        case 'Forecast Semicond'
            figmain = 20000;
    end
end

			
            
