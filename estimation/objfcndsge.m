function [obj, varargout] = objfcndsge(para,YY,YY0,nobs,nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,para_mask,...
  para_fix,marglh,coint,cointadd,cointall,YYcoint0,MIN, nant, antlags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% obj evaluates the log  posterior = log likelihood + log prior
%% MIN takes into account that csminwel minimizes as opposed to max, and that the parameters are rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MIN
    para = trans(para,trspec);
end
para = para.*(1-para_mask)+para_fix.*para_mask;


%--- Evaluate the Likelihood function: LOG{P[y|theta]}
lnpy = dsgelh(para,YY,YY0,nobs,nlags,nvar,mspec,npara,coint,cointadd,YYcoint0, nant, antlags);
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
