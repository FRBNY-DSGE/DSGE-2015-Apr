% OVERVIEW
%
% Function checks if a model is a two part model
function [ yesno ] = is2part(mspec)
  class2part;
  yesno = (any(mspec == class2part_all));
end
