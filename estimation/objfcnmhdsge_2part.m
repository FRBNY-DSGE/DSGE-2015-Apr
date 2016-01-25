function [lnpost,lnpy,zend,ZZ,DD,QQ] = objfcnmhdsge_2part(para,bounds,YY,...
    YY0,nobs,nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,...
    TTT,RRR,CCC,valid,para_mask,coint,cointadd,cointall,YYcoint0,nant,antlags,nshocks)

iSplit = getState(mspec,0,'n_end+n_exo');
nstate = size(TTT,1);

% For 510,555, TTT/RRR/CCC for 510 are contained within those for 555
if any(mspec == [555 556 557 5571 558])
    TTT1 = zeros(nstate-(nant+1),nstate-(nant+1));  
    TTT1 = [TTT(1:iSplit,1:iSplit), ...
             TTT(1:iSplit,iSplit+(nant+1)+1:end);...
             TTT(iSplit+(nant+1)+1:end,1:iSplit),...
             TTT(iSplit+(nant+1)+1:end,iSplit+(nant+1)+1:end)];
    
    RRR1 = zeros(nstate-(nant+1),nshocks);
    RRR1 = [RRR(1:iSplit,1:nshocks); RRR(iSplit+(nant+1)+1:end,1:nshocks)];
    
    CCC1 = zeros(nstate-(nant+1),1);
    CCC1 = [CCC(1:iSplit,1); CCC(iSplit+(nant+1)+1:end,1)];
end

if valid < 1;
    lnpost = -1E10;
    lnpy   = -1E20;
    zend = [];
    ZZ = [];
    DD = [];
    QQ = [];
    return;
end

ind1 = (para > bounds(:,2));
%disp([find(ind1) para(find(ind1)) bounds(find(ind1),1)])
ind1(logical(para_mask)) = [];
ind2 = (para < bounds(:,1));
%disp([find(ind2) disp(find(ind2)) bounds(find(ind2),2)])
ind2(logical(para_mask)) = [];

if any(ind1) || any(ind2)

    lnpost = -1E10;
    lnpy   = -1E20;
    zend = [];
    ZZ = [];
    DD = [];
    QQ = [];
    
else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: define the measurement equation: X_t = ZZ S_t +D+u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(strcat('[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT1,RRR1,valid,para,nvar-nant,nlags,mspec,npara,coint,cointadd);'));
%[ZZ,DD,QQ,EE,MM,retcode] = measur(TTT,RRR,valid,para,nvar,nlags,mspec,npara);

if retcode == 0
    % invalid parameterization
    pyt = -1E10;
    return;
    
end;

nstate  = size(TTT1,1);
nshocks = size(RRR1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: compute log-likelihood using Kalman filter - written by Iskander
%%         note that Iskander's program assumes a transition equation written as:
%%         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t 
%%         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
%%         and  VV2 = cov(eps2_t,u_u) = RRR*VV
%%         define VVall as the joint variance of the two shocks VVall = var([eps2_t;u_t])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

if coint == 0,

    % truncate data
    YY1 = YY(:,1:nvar-nant);

    HH = EE+MM*QQ*MM';
    VV = QQ*MM';
    VVall = [[RRR1*QQ*RRR1',RRR1*VV];[VV'*RRR1',HH]];

    if any(CCC ~= 0)
      DDadd = ZZ*((eye(size(TTT1))-TTT1)\CCC1);
    else
      DDadd = 0;
    end
    
    DD= DD+DDadd;

    %% if you have constant, change the st st of observables and rewrite kalman as deviations from stst


    %% Define the initial mean and variance for the state vector
    % (deals with nonstationarity on model by model basis)
    [A0,P0] = lyap_nonstationary(mspec,para,TTT1,RRR1,QQ);


    if ~isempty(YY0)
    YY0 = YY0(:,1:nvar-nant);
    %% run Kalman filter on initial observations
    [xxx,zend,Pend] = ...
        kalcvf2NaN(YY0',1,zeros(nstate,1),TTT1,DD,ZZ,VVall,A0,P0);

    else
      zend = A0;
      Pend = P0;
    end

    %% run Kalman filter on main sample  
    % First, filter using normal model up to period 
    [pyt1,zend,Pend] = kalcvf2NaN(YY1(1:end-antlags-1,:)',1,zeros(nstate,1),TTT1,DD,ZZ,VVall,zend,Pend);
    
    % Change models to incorporate anticipated shocks
    eval(strcat('[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT,RRR,valid,para,nvar,nlags,mspec,npara,coint,cointadd,nant);'));

    HH = EE+MM*QQ*MM';
    VV = QQ*MM';
    VVall = [[RRR*QQ*RRR',RRR*VV];[VV'*RRR',HH]];
    
    nstate  = size(TTT,1);
    nshocks = size(RRR,2);
    
    [zprev,pprev] = augmentStates(mspec,nant,nstate,zend,Pend);
    [pyt2,zend,pend]  = kalcvf2NaN(YY(end-antlags:end,:)',1,zeros(nstate,1),TTT,DD,ZZ,VVall,zprev,pprev);
    pyt=pyt1+pyt2;

else
    
    %% use system matrices for extended system
       
       if find(CCC) ~= nan,
           dfgdfsgdfgdf % if we get this error, we need to fix something
       end;
       
       DDcoint = DD(nvar+1:nvar+coint,:);
       ZZcoint = ZZ(nvar+1:nvar+coint,:);

       HH = EE+MM*QQ*MM';
       VV = QQ*MM';
       VVall = [[RRR*QQ*RRR',RRR*VV];[VV'*RRR',HH]];
       
       %% if you have constant, change the st st of observables and rewrite kalman as deviations from stst
      % DDadd = ZZ*((eye(size(TTT))-TTT)\CCC);
      % DD= DD+DDadd;

       %% Define the initial mean and variance for the state vector
       A0 = zeros(nstate,1);
       P0 = dlyap(TTT,RRR*QQ*RRR');

       %% run Kalman filter for initial observations on extended system
            [pytcoint,zend,Pend] = ...
               kalcvf2NaN([YY0(1,:),YYcoint0]',1,zeros(nstate,1),TTT,DD,ZZ,VVall,A0,P0);

       %% run Kalman filter on remainder of initialization sample                              
       %% first, remove cointegration component from the system
       ZZ = ZZ(1:nvar,:);
       DD = DD(1:nvar,:);
       EE = EE(1:nvar,1:nvar);
       MM = MM(1:nvar,:);

       HH = EE+MM*QQ*MM';
       VV = QQ*MM';
       VVall = [[RRR*QQ*RRR',RRR*VV];[VV'*RRR',HH]];
       nstate = size(TTT,1);

       if nlags > 1
          [pyt0,zend,Pend] = ...
               kalcvf2NaN(YY0(2:end,:)',1,zeros(nstate,1),TTT,DD,ZZ,VVall,zend,Pend);
       end                              

       %% run Kalman filter on main sample                              
       [pyt,zend,Pend] = ...
             kalcvf2NaN(YY',1,zeros(nstate,1),TTT,DD,ZZ,VVall,zend,Pend);
                                                                                           
end;   

    lnpy    = real(pyt);
    lnprio  = priodens(para,pmean,pstdd,pshape);
    lnpost  = lnpy + real(lnprio);
    
end
