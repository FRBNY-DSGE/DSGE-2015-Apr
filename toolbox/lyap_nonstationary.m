function [A0,P0] = lyap_nonstationary(mspec,para,T,R,Q)

switch mspec
    case {557 5571 558}
        i = getState(mspec,0,'zlev_t');
        A0 = zeros(size(T,1),1);
        A0(i) = para(5);
        
        T_tr = T([1:i-1 i+1:end],[1:i-1 i+1:end]);
        R_tr = R([1:i-1 i+1:end],:);
        
        P0_tr = dlyap(T_tr,R_tr*Q*R_tr');
        
        P0 = zeros(size(T));
        
        P0(1:i-1,1:i-1) = P0_tr(1:i-1,1:i-1);
        P0(1:i-1,i+1:end) = P0_tr(1:i-1,i:end);
        P0(i+1:end,1:i-1) = P0_tr(i:end,1:i-1);
        P0(i+1:end,i+1:end) = P0_tr(i:end,i:end);
    case {904, 9043, 90451}
        if mspec == 90451
            i1 = getState(mspec,0,'ztil_t');
            A0 = zeros(size(T,1),1);


            [nrT, ncT] = size(T);
            nrR = size(R,1);

            irT = 1:nrT;
            irT([i1]) = [];
            icT = 1:ncT;
            icT([i1]) = [];

            irR = 1:nrR;
            irR([i1]) = [];

            T_tr = T(irT, icT);
            R_tr = R(irR,:);
            
            P0_tr = dlyap(T_tr,R_tr*Q*R_tr');
            
            P0 = zeros(size(T));

            P0(irT, icT) = P0_tr;
        else 
            A0 = zeros(size(T,1),1);
            P0 = dlyap(T,R*Q*R');
        end
    otherwise
        A0 = zeros(size(T,1),1);
        P0 = dlyap(T,R*Q*R');
end
        
