function [peachdataimplied_s] = getpeachdataimplied_s(npeach,A,total_smth)

% This program isolates part of the forecast code that was enraging parfor.

peachdataimplied_s = zeros(npeach,size(A,1));

for ti = 1:npeach
    peachdataimplied_s(ti,:) = A*total_smth(:,ti);
end