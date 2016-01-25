% OVERVIEW
%
% This function expands the number of states to accomodate extra states for the
% anticipated policy shocks. It does so by taking the zend and Pend for the
% state space without anticipated policy shocks, then shoves in nant (or
% nant+1, see below) zeros in the middle of zend and Pend in the location of
% the anticipated shock entries.
%
% Note, 557 needs to add in a new equation and state for the monetary shock (bc
% it's iid without ant), while 904 doesn't. As a result, 557 needs nant+1 extra
% states/zeros, while 955 needs nant. See dsgelh_getNoZB.m for more detailed
% explanation.

function [zprev, pprev] = augmentStates(mspec,nant,nstate,zend,Pend)

class2part;

% Because money shock named differently in different models; also need to
% specify "adj" -- the number of extra equations beyond nant that are added to
% the system 
if any(mspec == [904 class2part_955])
    r_tl1 = getState(mspec,nant,'rm_tl1');
else
    r_tl1 = getState(mspec,nant,'r_tl1');
end
adj = adj2part(mspec);

up_to_ant_all = 1:(r_tl1-1-adj); % Indices for all states up to (but not including) ant and (maybe) new monetary shock states
aftr_ant_prev = (r_tl1+nant):nstate; % State indices for everything after ant in the new zprev vector 
aftr_ant_ends = (r_tl1-adj):length(zend); % State indices for everything after ant in the old zend vector

zprev = zeros(nstate,1);
zprev(up_to_ant_all) = zend(up_to_ant_all, end);
zprev(aftr_ant_prev) = zend(aftr_ant_ends, end); 

pprev = zeros([nstate nstate]);
pprev(up_to_ant_all, up_to_ant_all) = Pend(up_to_ant_all, up_to_ant_all, end); % Upper left block
pprev(aftr_ant_prev, up_to_ant_all) = Pend(aftr_ant_ends, up_to_ant_all, end); % Lower left block
pprev(up_to_ant_all, aftr_ant_prev) = Pend(up_to_ant_all, aftr_ant_ends, end); % Upper right block
pprev(aftr_ant_prev, aftr_ant_prev) = Pend(aftr_ant_ends, aftr_ant_ends, end); % Lower right block


