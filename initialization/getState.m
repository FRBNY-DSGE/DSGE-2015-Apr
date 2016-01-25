% OVERVIEW
% getState.m gets the numerical index for a state. This function replaces
% getE_pi.m, getpigap_t.m, getr_sh.m, getrt.m, getR_t.m
%
% INPUTS
% name: the state name as a string
%
% OUTPUTS
% stateNum: state's numerical index
%
% EXAMPLE: 
% r_sh = getState(mspec,nant,'r_sh')
% returns r_sh = 10

function stateNum = getState(mspec,nant,namer)
eval(['states' num2str(mspec)]);
eval(['stateNum = ',namer,';']);