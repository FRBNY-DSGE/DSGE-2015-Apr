function [counter_new,counter_states] = getcounter_s(ShockremoveList,peachflag,psize,Enddate_forecastfile,...
                                  sm_states,sm_shocks,sm_states_peach,sm_shocks_peach,...
                                  Startdate,nobs,nshocks,qahead,nvar,TTT,RRR,DD,ZZ,ind_r,ind_r_sh,bdd_int_rate,A,nplotstates,varargin)

                              
if nargin == 25
    mspec = varargin{1};
end
counter_new = zeros(nvar*(Enddate_forecastfile-Startdate),peachflag+1,length(ShockremoveList));
counter_states = zeros(nplotstates*(Enddate_forecastfile-Startdate),peachflag+1,length(ShockremoveList));

for peachcount = 0:peachflag

    Enddate1 = Enddate_forecastfile;

    switch peachcount
        case 0 % Unconditional
            z0 = sm_states(:,Startdate);
            shocks_all = sm_shocks(:,Startdate+1:nobs)';
        case 1 % Full conditional (conditional on T+1 data for almost all variables)
            z0 = sm_states_peach(:,Startdate);
            shocks_all = sm_shocks_peach(:,Startdate+1:nobs+psize)';
    end
    
    shockcount = 0;
    
    for Shockremove = ShockremoveList
        shockcount = shockcount + 1;
        
        z = z0;
        
        shocks = shocks_all;

        if any(Shockremove == [1:nshocks]),
            shocks(:,Shockremove) = 0;
        elseif Shockremove == 0
            shocks(:,:) = 0;
        end
        
        % shocks(Enddate1-Startdate-qahead_extra+1:Enddate1-Startdate,:) = 0;
        
        switch peachcount
            case 0
                shocks(Enddate1-Startdate-qahead+1:Enddate1-Startdate,:) = 0; % JC 7/9/2012
            case {1,2}
                shocks(Enddate1-Startdate-(qahead-psize)+1:Enddate1-Startdate,:) = 0; % JC 7/9/2012
        end
        
        yycounter = zeros(Enddate1-Startdate,nvar);
        counter_s = zeros(Enddate1-Startdate,nplotstates);
        for t = 1:Enddate1-Startdate
            z_test = TTT*z+RRR*shocks(t,:)';
            yycounter(t,:) = (DD+ZZ*z_test)';
            counter_s(t,:) = (A*z_test)';


            % This changes the monetary policy shock to account for the 0.25 interest rate bound

            
            if ~exist('mspec','var') || ~any(mspec==[803 805 904 828 825 826 827 829 830 835 8351 8352 836 8361 8362 906 8281:8285])
                ZeroBound=0.25;
                ZeroBoundTest=0.249;
            else
                ZeroBound=0.25/4;
                ZeroBoundTest=ZeroBound-0.01;
            end
            
            if bdd_int_rate && (yycounter(t,ind_r)<ZeroBound)
                shocks(t,ind_r_sh) = 0;
                shocks(t,ind_r_sh) = 1/(ZZ(ind_r,:)*RRR(:,ind_r_sh))*(ZeroBound - DD(ind_r) - ZZ(ind_r,:)*(TTT*z + RRR*shocks(t,:)'));
                z = TTT*z+RRR*shocks(t,:)';
                yycounter(t,:) = (DD+ZZ*z)';
                counter_s(t,:) = (A*z)';
                if yycounter(t,ind_r) < ZeroBoundTest
                    error('0.25 interest rate bound procedure not working')
                end
            else
                z = z_test;
            end
        end
        
        counter_new(:,peachcount+1,shockcount) = yycounter(:)';
        counter_states(:,peachcount+1,shockcount) = counter_s(:)';
        
    end
end