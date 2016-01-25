% hpdata.m
% Larry Christiano
% June 21, 1990

function [yf,yt] = hpfilter(y,lam);

% This program applies the hp filter to data

d = length(y);
a = sparse(d,d);
for i = 3:d-2
	a(i,i) = 6*lam+1;
	a(i,i+1) = -4*lam;
	a(i,i+2) = lam;
	a(i,i-2) = lam;
	a(i,i-1) = -4*lam;
	end

a(2,2) = 1+5*lam;
a(2,3) = -4*lam;
a(2,4) = lam;
a(2,1) = -2*lam;
a(1,1) = 1+lam;
a(1,2) = -2*lam;
a(1,3) = lam;

a(d-1,d-1) = 1+5*lam;
a(d-1,d-2) = -4*lam;
a(d-1,d-3) = lam;
a(d-1,d) = -2*lam;
a(d,d) = 1+lam;
a(d,d-1) = -2*lam;
a(d,d-2) = lam;

yt = a\y;
yf = y-yt;

