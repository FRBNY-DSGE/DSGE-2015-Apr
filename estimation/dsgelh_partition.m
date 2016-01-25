% OVERVIEW
% 
% This function will partition the sample for estimation purposes,
% accommodating 2-part estimation if relevant.
%
% Returns a structure array where size = number of distinct periods (presample,
% normal, ZB, etc.). The mt struture will hold the relevant matrices for each
% period.
function [mt, pd] = dsgelh_partition(YY0, YY, nvar, nant, antlags)

% Indicates 2-part estimation
if nant > 0
  if ~isempty(YY0), YY0 = YY0(:,1:nvar-nant); end % Take out expectations
  mt = struct('YY', {YY0, YY(1:end-antlags-1, 1:nvar-nant), YY(end-antlags:end,:)}, ...
              'nvar', {nvar-nant, nvar-nant, nvar}, ... % Effective number of variables in each period
              'nant', {{}, {}, {nant}}); % Number of anticipated shocks in each period
  pd = 3; % Set number of periods; here, pre-sample, normal, ZB
else
  mt = struct('YY',   {YY0, YY}, ...
              'nvar', {nvar, nvar}, ...
              'nant', {{}, {}});
  pd = 2; % Set number of periods; here, pre-sample, normal
end


end
