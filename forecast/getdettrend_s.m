function [dettrend_all] = getdettrend_s(Enddate_forecastfile,nvar,ZZ_,TTT_,DD,z0T,Startdate,peachflag,A,nplotstates)

dettrend_all = zeros((Enddate_forecastfile-Startdate+1)*nvar,peachflag+1);
dettrend_states = zeros((Enddate_forecastfile-Startdate+1)*nplotstates,peachflag+1);

for peachcount = 0:peachflag
    
    z = z0T(:,peachcount+1);
    
    dettrend = zeros(Enddate_forecastfile,nvar);
    dettrend_s = zeros(Enddate_forecastfile,nplotstates);
    for t = 1:Enddate_forecastfile
        dettrend(t,:) = (DD + ZZ_*(TTT_^t)*z)';
        dettrend_s(t,:) = (A*(TTT_^t)*z)';
    end

    dettrend_short = dettrend(Startdate:end,:);
    dettrend_s_short = dettrend_s(Startdate:end,:);

    dettrend_all(:,peachcount+1) = dettrend_short(:)';
    dettrend_states(:,peachcount+1) = dettrend_s_short(:)';
end
