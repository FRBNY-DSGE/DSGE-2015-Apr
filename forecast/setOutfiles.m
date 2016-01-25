% OVERVIEW
% setOutfiles.m creates variable outfile, where each field name is a type of
%               data that we want to write (ie forecast or shockdec).
%               Each field of outfile is the name of the output file.
%
% IMPORTANT VARIABLES
% outfile: a structure where each field name is a type of
%          data that we want to write. Each field of outfile is the
%          name of the output file. Example: outfile.forecast =
%          '/data/data_dsge_dir/.../forecastm...."
% outdataType: a cell array listing the type of data you want to create
%              outfiles for. Example: {forecast,datashocks}.

function [outfile] = setOutfiles(outdataType,zerobound, spath,peachflag,ShockremoveList)

for iDT = 1:length(outdataType)
    switch outdataType{iDT}
        case 'forecast'
          if zerobound
            outfile.forecast = [spath,'/forecast'];
            outfile.forecastbdd = [spath,'/forecastbdd'];
          else
            outfile.forecast = [spath,'/forecast_nozb'];
          end
        case 'condforecast'
          if zerobound
            outfile.condforecast = [spath,'/condforecast'];
            outfile.condforecastbdd = [spath,'/condforecastbdd'];
          else
            outfile.condforecast = [spath,'/condforecast_nozb'];
          end
        case 'states'
            outfile.states = [spath,'/states'];
        case 'peachshocks'
            outfile.peachshocks = [spath,'/peachshocks'];
        case 'semipeachshocks'
            outfile.semipeachshocks = [spath,'/semipeachshocks'];
        case 'condshocks'
            outfile.condshocks = [spath,'/condshocks'];
        case 'semicondshocks'
            outfile.semicondshocks = [spath,'/semicondshocks'];
        case 'datashocks'
            outfile.datashocks = [spath,'/datashocks'];
        case 'peachshocks_ns'
            outfile.peachshocks_ns = [spath,'/peachshocks_ns'];
        case 'semipeachshocks_ns'
            outfile.semipeachshocks_ns = [spath,'/semipeachshocks_ns'];
        case 'condshocks_ns'
            outfile.condshocks_ns = [spath,'/condshocks_ns'];
        case 'semicondshocks_ns'
            outfile.semicondshocks_ns = [spath,'/semicondshocks_ns'];
        case 'datashocks_ns'
            outfile.datashocks_ns = [spath,'/datashocks_ns'];
        case 'counter'
          for peachcounter = 0:peachflag
            for Shockremove = ShockremoveList
              switch peachcounter
                case 1, outfile.(['counter',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/counter_peach', num2str(Shockremove)];
                case 0, outfile.(['counter',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/counter', num2str(Shockremove)];
              end
            end
          end
        case 'shockdec'
          for peachcounter = 0:peachflag
            for Shockremove = ShockremoveList
              switch peachcounter
                case 1, if Shockremove > 0, outfile.(['shockdec',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/shockdec_peach', num2str(Shockremove)]; end
                case 0, if Shockremove > 0, outfile.(['shockdec',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/shockdec', num2str(Shockremove)]; end
              end
            end
          end
        case 'ytrend'
            outfile.ytrend = [spath,'/ytrend'];
        case 'dettrend'
            outfile.dettrend = [spath,'/dettrend'];
        case 'dettrend_peach'
            outfile.dettrend_peach = [spath,'/dettrend_peach'];
        case 'hist'
            outfile.hist = [spath,'/hist'];
        case 'condhist'
            outfile.condhist = [spath,'/condhist'];
        case 'forecastPRIO'
            outfile.forecast = [spath,'/forecastPRIO'];
            
        case 'forecast_s'
          if zerobound
            outfile.forecast_s = [spath,'/forecast_s'];
            outfile.forecastbdd_s = [spath,'/forecastbdd_s'];
          else
            outfile.forecast_s = [spath,'/forecast_s'];
          end
        case 'condforecast_s'
          if zerobound
              outfile.condforecast_s = [spath,'/condforecast_s'];
              outfile.condforecastbdd_s = [spath,'/condforecastbdd_s'];
          else
              outfile.condforecast_s = [spath,'/condforecast_s'];
          end
            
        case 'dettrend_s'
            outfile.dettrend_s = [spath,'/dettrend_s'];
        case 'dettrend_peach_s'
            outfile.dettrend_peach_s = [spath,'/dettrend_peach_s'];
        case 'shockdec_s'
          for peachcounter = 0:peachflag
            for Shockremove = ShockremoveList
              switch peachcounter
                case 1, if Shockremove > 0, outfile.(['shockdec_s',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/shockdec_peach_s', num2str(Shockremove)]; end
                case 0, if Shockremove > 0, outfile.(['shockdec_s',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/shockdec_s', num2str(Shockremove)]; end
              end
            end
          end
        case 'counter_s'
          for peachcounter = 0:peachflag
            for Shockremove = ShockremoveList
              switch peachcounter
                case 1, outfile.(['counter_s',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/counter_peach_s', num2str(Shockremove)];
                case 0, outfile.(['counter_s',num2str(peachcounter),num2str(Shockremove)]) = [spath,'/counter_s', num2str(Shockremove)];
              end
            end
          end
        case 'hist_s'
            outfile.hist_s = [spath,'/hist_s'];
        case 'hist_peach_s'
            outfile.hist_peach_s = [spath,'/hist_peach_s'];
    end
end
