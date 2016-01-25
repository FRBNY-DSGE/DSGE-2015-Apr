function Ys = getYs(YY_allvars,YY_p_allvars,dlpop,dlpop_p,cum_for,popadj,nvar,fourq_flag,q_adj,noadj,varargin)

if nargin == 11
    pop_smth = varargin{1};
end

Ys = zeros(size(YY_allvars));

if noadj
    YY = YY_allvars(:,I);
    YY_p = YY_p_allvars(:,I);

    if cum_for(I) == 2
        YYold = [YY_p(end); YY(1:end-1)];
        Ys(:,I) = YY- YYold;
    else
        Ys(:,I) = YY;
    end

else

    % Commented out bc annoying AF
    %if q_adj ~= 400, warning('q_adj is not 400'); end

    % May need to agument YY_p with nans so that the four quarter
    % calculations go through
    N_p = size(YY_p_allvars,1);
    if fourq_flag && N_p < 4
      YY_p_allvars = [nan(4-N_p, size(YY_p_allvars,2)); YY_p_allvars];
      dlpop_p = [nan(4-N_p, size(dlpop_p,2)); dlpop_p];
    end


    for I = 1:nvar

        YY = YY_allvars(:,I);
        YY_p = YY_p_allvars(:,I);

        if cum_for(I) == 1 && popadj(I) == 0

            if fourq_flag
                % Transform from log growth rates to 4Q growth rates.

                YY1 = [YY_p(end); YY(1:end-1)];
                YY2 = [YY_p(end-1:end); YY(1:end-2)];
                YY3 = [YY_p(end-2:end); YY(1:end-3)];

                Ys(:,I) = 100*(exp((YY + YY1+ YY2 + YY3)/q_adj) - 1);

            else
                % Transform from log growth rates to quarterly annualized
                % growth rates.

                Ys(:,I) = 100*(exp(YY/q_adj).^4 - 1);

            end

        elseif cum_for(I) == 1 && popadj(I) == 1

            if fourq_flag
                % Transform from log per capita growth rates to 4Q total
                % growth rates.

                YYpop = YY/q_adj + dlpop;
                YYpop1 = [YY_p(end); YY(1:end-1)]/q_adj + [dlpop_p(end); dlpop(1:end-1)];
                YYpop2 = [YY_p(end-1:end); YY(1:end-2)]/q_adj + [dlpop_p(end-1:end); dlpop(1:end-2)];
                YYpop3 = [YY_p(end-2:end); YY(1:end-3)]/q_adj + [dlpop_p(end-2:end); dlpop(1:end-3)];

                Ys(:,I) = 100*(exp(YYpop + YYpop1 + YYpop2 + YYpop3) - 1);

            else
                % Transform from log per capita growth rates to quarterly
                % annualized total growth rates.

                Ys(:,I) = 100*((exp(YY/q_adj + dlpop)).^4-1);

            end

        elseif cum_for(I) == 2 && popadj(I) == 0

            if fourq_flag
                % Transform from log levels to 4Q growth rates.

                YYold = [YY_p(end-3:end); YY(1:end-4)];
                Ys(:,I) = 100*(exp(YY/100 - YYold/100) - 1);

            else
                % Transform from log levels to quarterly annualized growth
                % rates.

                YYold = [YY_p(end); YY(1:end-1)];
                Ys(:,I) = 100*((exp(YY/100 - YYold/100).^4) - 1);

            end

        elseif cum_for(I) == 2 && popadj(I) == 1

            if fourq_flag
                % Transform from log per capita levels to 4Q total growth
                % rates.

                dlpop1 = [dlpop_p(end); dlpop(1:end-1)];
                dlpop2 = [dlpop_p(end-1:end); dlpop(1:end-2)];
                dlpop3 = [dlpop_p(end-2:end); dlpop(1:end-3)];

                dlpop_4q = dlpop + dlpop1 + dlpop2 + dlpop3;

                YYold = [YY_p(end-3:end); YY(1:end-4)];

                Ys(:,I) = 100*(exp(YY/100 - YYold/100 + dlpop_4q) - 1);

            else
                % Transform from log per capita levels to quarterly annualized
                % total growth rates.

                YYold = [YY_p(end); YY(1:end-1)];

                Ys(:,I) = 100*(exp(YY/100 - YYold/100 + dlpop).^4 - 1);

            end

        elseif cum_for(I) == 3
            % Transform from log levels to levels

            Ys(:,I) = exp(YY/100);

        elseif cum_for(I) == 4
            %This is for unemployment
            Ys(:,I)=(exp(YY/100)-1)*100;
        elseif cum_for(I) == 5
            Ys(:,I) = YY*4;
        elseif cum_for(I) ==6 && popadj(I) ==1
            Ys(:,I) = YY + dlpop*100;
        elseif cum_for(I) ==6 && popadj(I) ==0
            Ys(:,I) = YY;

        elseif cum_for(I) == 7
            % Output is in log levels
            Ys(:,I) = YY;
        elseif cum_for(I) ==0 && popadj(I) ==1
            Ys(:,I) = YY + 100*log(pop_smth);
        else

            Ys(:,I) = YY*400/q_adj;

        end

    end

end