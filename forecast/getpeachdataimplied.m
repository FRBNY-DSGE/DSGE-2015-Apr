function [peachdataimplied] = getpeachdataimplied(npeach,ZZ_,DD,total_smth)

% This program isolates part of the forecast code that was enraging parfor.

peachdataimplied = zeros(npeach,size(ZZ_,1));

for ti = 1:npeach
    peachdataimplied(ti,:) = (ZZ_*total_smth(:,ti) + DD)';
end