function [Vseq,Vseq_alt,Shseq,Shseq_alt] = setVseq(plotType,nshocks,nant,nvar,mspec,varargin)

class2part
adj = 0;

%% Preallocate space for memory usage
Vseq=[];
Vseq_alt=[];
Shseq=[];
Shseq_alt=[];

%% Settings
  switch plotType
      case {'Forecast','Shock Decomposition','Shock Decomposition-noDet','Exp_Forecast'}
            if any(mspec == class2part_all)
              Vseq = 1:(nvar-nant+adj);
              Vseq_alt = [];
            else 
              Vseq = 1:nvar;
              Vseq_alt = [];
            end
      case {'Forecast Semicond'}
          Vseq = [1 5 6];
          Vseq_alt = [];
          Vseq = [Vseq Vseq_alt];
      case {'Forward Guidance'}
          Vseq = [14 15 16 17 12 13];
          Vseq_alt = [];
          Vseq = [Vseq Vseq_alt];
      case 'Forward Guidance Horizontal'
          Vseq = [14 15 12 16 17 13];
          Vseq_alt = [];
          Vseq = [Vseq Vseq_alt];
      case {'Counterfactual by Variable'}
              Shseq = 1:nshocks-nant;
              Shseq_alt = [];
              if any(mspec == class2part_all)
                Vseq = 1:(nvar-nant+adj); 
              else 
                Vseq = 1:nvar;
              end
              Vseq_alt = [];
      case {'Counterfactual by Shock'}
              if any(mspec == class2part_all)
                Shseq = 1:(nvar-nant+adj);
              else 
                Shseq = 1:nvar;
              end
              Shseq_alt = [];
              Vseq = 1:nshocks-nant; % These are now shocks
              Vseq_alt = [];

      case {'Shock','Shock Squared','Htil','Eta Squared','Sigma'}
            if any(mspec == class2part_all)
                Vseq = 1:nshocks-nant; % These are now shocks
            else
                Vseq = 1:nshocks; % These are now shocks
            end
            Vseq_alt = [];

      case {'ShockAnt'}
              if is2part(mspec)
                  Vseq = [nshocks-nant+1:nshocks];
                  Vseq_alt = [];
              else
                  Vseq = 1:nshocks;
                  Vseq_alt = [];
              end

      case {'Forecast Comparison', '4Q Forecast', 'Forecast Comparison Separate', '4Q Forecast Separate'}

              Vseq = 1:nvar;
              Vseq_alt = [];

          
      case 'ImpObs'
          nimpVar = varargin{1};
          Vseq = nvar+1:nvar+nimpVar;
  end
  
  