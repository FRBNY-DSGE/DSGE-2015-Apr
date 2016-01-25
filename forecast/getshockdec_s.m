function [shockdec_new,shockdec_states] = getshockdec_s(ShockremoveList,peachflag,psize,Enddate_forecastfile,...
                                    sm_shocks,sm_shocks_peach,...
                                    Startdate,nobs,qahead,nvar,TTT,RRR,DD,ZZ,ind_r,ind_r_sh,bdd_int_rate,A,nplotstates)

shockdec_new = zeros(nvar*(Enddate_forecastfile-Startdate+1),peachflag+1,length(ShockremoveList)-1);
shockdec_states = zeros(nplotstates*(Enddate_forecastfile-Startdate+1),peachflag+1,length(ShockremoveList)-1);

for peachcount = 0:peachflag

    Enddate1 = Enddate_forecastfile;
    
    switch peachcount
        case 0 % Unconditional
            shocks_all = sm_shocks(:,1:nobs)';
        case 1 % Full conditional (conditional on T+1 data for almost all variables)
            shocks_all = sm_shocks_peach(:,1:nobs + psize)';

    end
    
    if peachcount > 0

    else

    end

    shockcount = 0;
    
    for Shockremove = ShockremoveList
        
        if Shockremove > 0
            
            shockcount = shockcount + 1;
            
            z = zeros(length(TTT),1);

            shocks = zeros(size(shocks_all));
            shocks(:,Shockremove) = shocks_all(:,Shockremove);

            switch peachcount
                case 0
                    shocks(Enddate1 - (qahead) + 1:Enddate1,:) = 0; % JC 7/9/2012
                case {1 2}
                    shocks(Enddate1 - (qahead - psize) + 1:Enddate1,:) = 0; % JC 7/9/2012
            end

            % Use RRR and not RRR1 as the shocks are non standardized
            yyshockdec = zeros(Enddate1,nvar);
            shockdec_s = zeros(Enddate1,nplotstates);
            for t = 1:Enddate1
                
                z = TTT*z + RRR*shocks(t,:)';
                yyshockdec(t,:) = (DD + ZZ*z)';
                shockdec_s(t,:) = (A*z)';
            end
            
            yyshockdec_short = yyshockdec(Startdate:Enddate1,:);
            shockdec_s_short = shockdec_s(Startdate:Enddate1,:);

            shockdec_new(:,peachcount+1,shockcount) = yyshockdec_short(:)';
            shockdec_states(:,peachcount+1,shockcount) = shockdec_s_short(:)';
        end
        
    end
    
end