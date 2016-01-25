function hist_s = gethistory_s(A,nplotstates,peachflag,stime,sm_states,sm_states_peach)

hist_s = zeros(nplotstates*stime,peachflag+1);

for peachcount = 0:peachflag
    
    if peachcount == 0
        states = sm_states(:,1:stime);
    elseif peachcount == 1
        states = sm_states_peach(:,1:stime);
    end
    tmp = zeros(stime,nplotstates);
    
    for I = 1:stime
        tmp(I,:) = (A*states(:,I))';
    end
    
    hist_s(:,peachcount+1) = tmp(:);
    
end