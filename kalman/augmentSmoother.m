% OVERVIEW
% augmentSmoother.m augments output matrices to accomodate 
%                   additional states for anticipated shocks. 
% 
% The state transition equation is: S_t = TTT*S_t-1+RRR*eps_t
% The measurement equation is: Y_t = ZZ*S_t+DD
%
% INPUTS
% r_tlx, r_tl1: 
% nant: number of anticipated shocks
% antlags: number of past periods in which the zerobound was in effect
% A0,P0: initial state vector and state covariance matrix
% pred0,pred,pred_ant: filtered state vector over (1) presample, (2) main sample up to
%                      T-antlags-1, and (3) main sample from T-antlags:T
% vpred0,vpred,vpred_ant: filtered state covariance matrix over (1) presample, 
%                         (2) main sample up to T-antlags-1, and (3) main sample from T-antlags:T
% YY0,YY,YY_ant: (1) pre-sample data (2) main sample data up to
%                T-antlags-1, and (3) main sample data from T-antlags:T

% OUTPUTS
% NOTE: *_all variables concatenate pre-sample and main sample including 
% anticipated shocks period
% pred_all: time series of the state vector 
% vpred_all: time series of the state covariance matrix 
% YY_all: time series of observable data
% A0_,P0_: augmented initial state vector and state covariance matrix

function [pred_all,vpred_all,YY_all,A0_,P0_] = augmentSmoother(nstate,r_tlx,r_tl1,nant,antlags,...
    pred0,pred,pred_ant,...
    vpred0,vpred,vpred_ant,...
    YY0,YY,YY_ant,...
    A0,P0)


pred_all = zeros(nstate,size(YY,1)+size(YY0,1));

if size(YY0,1) ~= 0
    pred_all(1:r_tlx-2,1:size(YY0,1)) = pred0(1:r_tlx-2,:);
    pred_all(r_tl1+nant:end,1:size(YY0,1)) = pred0(r_tlx-1:end,:);
end
pred_all(1:r_tlx-2,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = pred(1:r_tlx-2,:);
pred_all(r_tl1+nant:end,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = pred(r_tlx-1:end,:);

pred_all(:,size(YY0,1)+size(YY,1)-antlags:end) = pred_ant;

vpred_all = zeros(nstate,nstate,size(YY,1)+size(YY0,1));
if size(YY0,1) ~= 0 % this is because 557 doesn't work with presample
    vpred_all(1:r_tlx-2,1:r_tlx-2,1:size(YY0,1)) = vpred0(1:r_tlx-2,1:r_tlx-2,:);
    vpred_all(r_tl1+nant:end,1:r_tlx-2,1:size(YY0,1)) = vpred0(r_tlx-1:end,1:r_tlx-2,:);
    vpred_all(1:r_tlx-2,r_tl1+nant:end,1:size(YY0,1)) = vpred0(1:r_tlx-2,r_tlx-1:end,:);
    vpred_all(r_tl1+nant:end,r_tl1+nant:end,1:size(YY0,1)) = vpred0(r_tlx-1:end,r_tlx-1:end,:);
end
vpred_all(1:r_tlx-2,1:r_tlx-2,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = vpred(1:r_tlx-2,1:r_tlx-2,:);
vpred_all(r_tl1+nant:end,1:r_tlx-2,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = vpred(r_tlx-1:end,1:r_tlx-2,:);
vpred_all(1:r_tlx-2,r_tl1+nant:end,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = vpred(1:r_tlx-2,r_tlx-1:end,:);
vpred_all(r_tl1+nant:end,r_tl1+nant:end,size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1) = vpred(r_tlx-1:end,r_tlx-1:end,:);

vpred_all(:,:,size(YY,1)+size(YY0,1)-antlags:end) = vpred_ant;

% The data also needs to be concatated, with NaNs in
% the expectation terms for all periods prior to the
% model switch.

YY_all = NaN(size(YY,1)+size(YY0,1),size(YY_ant,2));

YY_all(1:size(YY0,1),1:size(YY0,2)) = YY0;
YY_all(size(YY0,1)+1:size(YY0,1)+size(YY,1)-antlags-1,1:size(YY,2)) = YY(1:end-antlags-1,:);
YY_all(end-antlags:end,:) = YY_ant;

% The starting expectation and covariance also need to
% be augmented (with zeros).

A0_ = zeros(nstate,1);
A0_(1:r_tlx-2) = A0(1:r_tlx-2);
A0_(r_tl1+nant:end) = A0(r_tlx-1:end);

P0_ = zeros(nstate);
P0_(1:r_tlx-2,1:r_tlx-2) = P0(1:r_tlx-2,1:r_tlx-2);
P0_(1:r_tlx-2,r_tl1+nant:end) = P0(1:r_tlx-2,r_tlx-1:end);
P0_(r_tl1+nant:end,1:r_tlx-2) = P0(r_tlx-1:end,1:r_tlx-2);
P0_(r_tl1+nant:end,r_tl1+nant:end) = P0(r_tlx-1:end,r_tlx-1:end);