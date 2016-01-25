function hpdband = hpdint(draws,percent,short);
%%  input: 
%%  draws: [ndraws x drawdim] matrix
%%  percent: peremtage WITHIN bands (say .9 to have 90% of mass within the bands) 
%%  short = 1, if choose shortest interval, otherwise just chop off lowest and highest percent/2
%%  output:
%%  hpbands: [2 x drawdim] matrix, first row: lower band, second row: upper band

[ndraws,drawdim]  = size(draws);
hpdband   = zeros(2,drawdim);
nwidth    = round(percent*ndraws);

for i = 1:drawdim;
	drawcoli = draws(:,i);
	% sort response for period i, element 1 is max
	drawcoli = flipud(sort(drawcoli));
    if short
		bup   = 1;
		minwidth  = drawcoli(1) - drawcoli(nwidth);
        done = 0;
        j = 2;
        while j <= (ndraws-nwidth+1)
        	newwidth = drawcoli(j) - drawcoli(j+nwidth-1);
		    if newwidth < minwidth;
                bup = j;
                minwidth = newwidth;                
            end;
        	j = j+1;
        end;

    else
        bup = ndraws-nwidth-floor(.5*(ndraws-nwidth));
    end

	hpdband(2,i) = drawcoli(bup);
	hpdband(1,i) = drawcoli(bup+nwidth-1);
end;
