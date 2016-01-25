% OVERVIEW
%
% Returns "adj", which is used to adjust the number of variables and sizing for
% 2part models. We need to do this because 555, 556, 557 add another equation
% and state to model a non-iid evolution of the monetary shock. On the other
% hand, 955 already has that
function [ adj ] = adj2part(mspec)
  class2part;
  if any(mspec == class2part_555)
    adj = 1;
  else
    adj = 0;
  end
end
