function dm(B)
%% fid = 1: print to screen
%% digit: after comma
fid = 1;
digit = 4;

format = '\n ';
for i = 1:size(B,2)
	format = strcat(format,' %2.',num2str(digit),'f ');
end

fprintf(fid,format,B');
fprintf(fid,'\n\n');

