function [Yaxis] = setYaxis(mspec,dataset,zerobound,nvar,plotType,V_a)

switch plotType
  case {'Forecast','Forecast Comparison','Exp_Forecast','Forward Guidance','Forward Guidance Horizontal','Forecast Semicond'}
    switch V_a
      case 1 	%% output
        Yaxis.lim = [-9 9];
        Yaxis.freq = 3;
      case 2	%% hours
        Yaxis.lim = [-9 12];
        Yaxis.freq = 1;
      case 3	%% labor share
          Yaxis.lim = [0.5 0.6];
          Yaxis.freq = 0.01;
      case 4	%% inflation
        Yaxis.lim = [0 4];
        Yaxis.freq = 1;
      case 5	%% interest rate
        if zerobound
          Yaxis.lim = [0 6];
          Yaxis.freq = 1;
        else
          Yaxis.lim = [-4 6];
          Yaxis.freq = 1;
        end
      case 6	
        Yaxis.lim = [-3 6];
        Yaxis.freq = 1;
      case 7	
        Yaxis.lim = [-10 10];
        Yaxis.freq = 2;
      case 8
        if zerobound
            Yaxis.lim = [-4 2];
            Yaxis.freq = 1;
        else
            Yaxis.lim = [-6 2];
            Yaxis.freq = 2;
        end
      case 9
        Yaxis.lim = [-8 8];
        Yaxis.freq = 2;
      case 10 
        Yaxis.lim = [-3 18];
        Yaxis.freq = 3;
      case 11
        Yaxis.lim = [-4 23];
        Yaxis.freq = 3;
      case 12 
        Yaxis.lim = [0 11];
        Yaxis.freq = 1;
      case 13	%% interest rate
        if zerobound
          Yaxis.lim = [0 6];
          Yaxis.freq = 1;
        else
          Yaxis.lim = [-4 6];
          Yaxis.freq = 1;
        end
      case 14 %% 4q output growth
          Yaxis.lim = [-5 13];
          Yaxis.freq = 1;
      case 15 %% Cumulative output growth
          Yaxis.lim = [-4 26];
          Yaxis.freq = 2;
      case 16
          Yaxis.lim = [0 5];
          Yaxis.freq = 0.5;
      case 17
          Yaxis.lim = [0 30];
          Yaxis.freq = 2;

    end
case {'4Q Forecast'}
    switch V_a
      case 1, 
        Yaxis.lim = [-6 10];
        Yaxis.freq = 1;
      case 2, 
        Yaxis.lim = [-15 6];
        Yaxis.freq = 4;
      case 3	%% labor share
        Yaxis.lim = [0.50 0.6];
        Yaxis.freq = 0.02;
      case 4
        Yaxis.lim = [-1 3.5];
        Yaxis.freq = 0.5;
      case 5
        if zerobound
            Yaxis.lim = [0 5.5];
            Yaxis.freq = 0.5;
        end
      case 6	
        Yaxis.lim = [-3 6];
        Yaxis.freq = 1;
      case 7	%% investment in mspec 6
        Yaxis.lim = [0 14];
        Yaxis.freq = 2;
      case 8
        if zerobound
          Yaxis.lim = [-4 2];
          Yaxis.freq = 1;
        else
          Yaxis.lim = [-6 2];
          Yaxis.freq = 2;
        end
      case 9
        Yaxis.lim = [-8 8];
        Yaxis.freq = 2;
      case 10 
        Yaxis.lim = [-3 18];
        Yaxis.freq = 3;
      case 11
        Yaxis.lim = [-4 23];
        Yaxis.freq = 3;
      case 12
          Yaxis.lim = [0 11];
          Yaxis.freq = 1;
      case 13	%% interest rate
        if zerobound
          Yaxis.lim = [0 6];
          Yaxis.freq = 1;
        else
          Yaxis.lim = [-4 6];
          Yaxis.freq = 1;
        end
      case 14 %% 4q output growth
          Yaxis.lim = [-5 13];
          Yaxis.freq = 1;
      case 15 %% Cumulative output growth
          Yaxis.lim = [-4 26];
          Yaxis.freq = 2;
      case 16
          Yaxis.lim = [0 5];
          Yaxis.freq = 0.5;
      case 17
          Yaxis.lim = [0 30];
          Yaxis.freq = 2;
    end
case {'Counterfactual','Counterfactual by Variable','Counterfactual by Shock'}
    switch V_a
        case 1 	%% output
            Yaxis.lim = [-10 10];
            Yaxis.freq = 4;
        case 2	%% hours
            Yaxis.lim = [-15 6];
%Yaxis.lim = [530 560]; % Use for levels
            Yaxis.freq = 4;
            %Yaxis.freq = 5; % Use for levels
        case 3	%% labor share
            Yaxis.lim = [0.50 0.6];
            Yaxis.freq = 0.02;
        case 4	%% inflation
            Yaxis.lim = [-1 5];
            Yaxis.freq = 2;
        case 5	%% interest rate
            Yaxis.lim = [-2 7];
            Yaxis.freq = 2;
        case 6	
            Yaxis.lim = [-3 6];
            Yaxis.freq = 1;
        case 7	
            Yaxis.lim = [-40 30];
            Yaxis.freq = 10;
        otherwise 
            Yaxis.lim = [-40 30];
            Yaxis.freq = 10;
    end
    
case {'Shock Decomposition','Shock Decomposition-noDet'}
  switch V_a
    case 1 	%% output
      Yaxis.freq = 2;
      Yaxis.lim = [-12 5];
    case 2	%% hours
      Yaxis.lim = [-9 6]; 
      Yaxis.freq = 1;
    case 3	%% labor share
      Yaxis.lim = [-4 3];
      Yaxis.freq = 0.5;
    case 4	%% inflation
      Yaxis.lim = [-2 2];
      Yaxis.freq = 1;
    case 5	%% interest rate
      Yaxis.lim = [-5 2];
      Yaxis.freq = 1;
    case 6
      Yaxis.lim = [-1 5];
      Yaxis.freq = 1;
    case 7
      Yaxis.lim = [-40 15];
      Yaxis.freq = 5;
    case 8
      Yaxis.lim = [-1 1];
      Yaxis.freq = 0.25;
  end
case {'Shock'}
  Yaxis.lim = [-4,4];
  if V_a==7, Yaxis.lim = [-7,7]; end
  Yaxis.freq = 0.5;
case {'ShockAnt'}
  Yaxis.lim = [-0.4,0.2];
  Yaxis.freq = 0.1;
case {'Shock Squared','Htil','Eta Squared','Sigma'}
  Yaxis.lim = [];
otherwise
  Yaxis.lim = [];
end

end
