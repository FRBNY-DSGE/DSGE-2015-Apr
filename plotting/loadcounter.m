function ypath_counter = loadcounter(fid_ctr,num_counter,nvar,counter_ahead,nsim,jstep)

ypath_counter = [];
while ( ftell(fid_ctr) < num_counter )
    % Read blocks of size nsim
    ypathadd_counter = fread(fid_ctr,[nvar*counter_ahead,nsim/jstep],'single')';
    ypath_counter = [ypath_counter;ypathadd_counter];
    ftell(fid_ctr);
end;
