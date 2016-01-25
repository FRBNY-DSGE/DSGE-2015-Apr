function [obj, varargout] = objfcndsge_2part(para,YY,nobs,nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,para_mask,...
  para_fix,marglh,coint,cointadd,cointall,MIN, nant, antlags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% obj evaluates the log  posterior = log likelihood + log prior
%% MIN takes into account that csminwel minimizes as opposed to max, and that the parameters are rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MIN
    para = trans(para,trspec);
end
para = para.*(1-para_mask)+para_fix.*para_mask;

%--- Evaluate the Likelihood function: LOG{P[y|theta]}
lnpy = dsgelh_2part(para,YY,nobs,nlags,nvar,mspec,npara,coint,cointadd,nant, antlags);
%--- Evaluate the Prior distribution: LOG{P[theta]}
lnprio = priodens(para,pmean,pstdd,pshape);

%--- Evaluate the Posterior density: LOG{P[theta|y]}
if MIN 
    obj = real(-lnpy-lnprio);  % We minize the inverse of the likelihood fcn 
else 
    obj = real(lnpy+lnprio);  
end

if nargout>1
    varargout(1) = {lnpy};
    varargout(2) = {lnprio};
end;
