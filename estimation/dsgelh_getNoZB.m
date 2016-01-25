% OVERVIEW
%
% This function pulls out the state equation matrices for the models WITHOUT
% anticipated policy shocks (model 510, 904, etc), from the state equation
% matrics for the models WITH anticipated policy shocks (e.g. 555, 557, 955).
% This is because the models WITH anticipated shocks contain the models WITHOUT
% anticipated shocks within them (i.e. 510 is embedded in 555 and 904 in 955).
function [ TTT, RRR, CCC ]  = dsgelh_getNoZB(nant, start_ant_state, start_ant_shock, TTT, RRR, CCC, varargin)

  if nargin > 6
    revol_ind = varargin{1};
    % This indicates that we need to drop an extra equation for the evolution
    % of the monetary shock. 
    %
    % We need to do this  because 555, 557 and variants don't normally have an
    % equation for the evolution of the contemporaneous monetary policy shock;
    % it's assumed iid and is implemented via PSI in the taylor rule equation
    % row. But in the model with anticipated policy shocks, an extra equation
    % has to be added (and a new state introduced) such that 
    % 
    %   (time t monetary shock) = (contemp shock) + (ant shocks). 
    %
    % This sum is then what's added into the Taylor rule. BUT if we kick out
    % the ant shocks, we don't need this equation, so let's remove it. We will
    % just treat the monetary shock as iid as in model 510.

%    TTT(:,revol_ind) = [];
%    TTT(revol_ind,:) = [];
%    RRR(revol_ind,:) = [];
%    CCC(revol_ind,:) = [];
  else
    revol_ind = [];

  end


  % Indices corresopnding to ant shock states (which hold exogenous ant shocks)
  % and exogenous shocks
  ant_state_inds = [revol_ind start_ant_state:(start_ant_state+nant-1)];
  ant_shock_inds = start_ant_shock:(start_ant_shock+nant-1);
  
  % Remove the ant rows and columns
  TTT(:,ant_state_inds) = []; % Ditch the columns that multiply the states that hold the anticipated policy shocks
  TTT(ant_state_inds,:) = []; % Ditch the equations for the evolution fo the anticipated policy shocks
  RRR(ant_state_inds,:) = []; 
  CCC(ant_state_inds,:) = []; 
  RRR(:,ant_shock_inds) = []; % Ditch columns for exogenous anticipated policy shocks


end
