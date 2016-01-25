% OVERVIEW
%
% Returns the starting indices for the anticipated policy shock states and
% shocks. Very model-class specific
function [start_ant_state, start_ant_shock, revol_ind] = get_start_ant(mspec, nant)
  
  class2part;

  % Set the starting indices for states and shocks corresponding to anticipated
  % policy shocks; depends on how eqconds is set up for a particular model 
  if any(mspec == class2part_555)
    str_ant_state = 'n_end+n_exo';
    str_ant_shock = 'n_exo';

    % Also for the 555 class of models, we need to drop the extra equation for
    % the evolution of the monetary shock. See the getNoZBmats.m function. This
    % is very model specific.
    revol_ind = {getState(mspec, nant, 'n_end+n_exo') - nant};

  elseif any(mspec == class2part_955)
    str_ant_state = 'nstates';
    str_ant_shock = 'nex';
    revol_ind = {};
  end
  start_ant_state = getState(mspec,nant,str_ant_state)-nant+1;
  start_ant_shock = getState(mspec,nant,str_ant_shock)-nant+1;
  
end
