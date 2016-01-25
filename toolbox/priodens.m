function lnprior = priodens(para,pmean,pstdd,pshape);
% This procedure computes a prior density for
% the structural parameters of the DSGE models
% pshape: 0 is point mass, both para and pstdd are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu)

  lnprior = 0;
  a = 0;
  b = 0;

  nprio = size(pshape,1);
  prioinfo = [zeros(nprio,2),pshape];

  for i = 1:nprio;
    if prioinfo(i,3) == 1; % BETA Prior %%check this distr
      a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
      b = a*(1/pmean(i) - 1);
      lnprior = lnprior + (a-1)*log(para(i))+(b-1)*log((1-para(i)))-betaln(a,b);   

    elseif prioinfo(i,3) == 2; % GAMMA PRIOR 
      b = pstdd(i)^2/pmean(i);
      a = pmean(i)/b;
      lnprior = lnprior + (a-1)*log(para(i))-(para(i)/b)-gammaln(a)-a*log(b);
    elseif prioinfo(i,3) == 3; % GAUSSIAN PRIOR 
      a = pmean(i);
      b = pstdd(i);
      lnprior = lnprior -0.5*log(2*pi)-log(b)-0.5*(para(i)-a)^2/b^2;
    elseif prioinfo(i,3) == 4; % INVGAMMA PRIOR 
      a = pmean(i);
      b = pstdd(i);
      %a= (.0264/b).^.5
      %para(i)
      lnprior = lnprior +log(2)-gammaln(b/2)+(b/2)*log(b*a^2/2)-((b+1)/2)*log(para(i)^2)-b*a^2/(2*para(i)^2);
      %v =log(2)-gammaln(b/2)+(b/2)*log(b*a^2/2)-((b+1)/2)*log(para(i)^2)-b*a^2/(2*para(i)^2)
    elseif prioinfo(i,3) == 5; % UNIFORM PRIOR pmean(i)=leftbound, pstdd(i)=rightbound
      a = pmean(i);% - sqrt(3)*pstdd(i);
      b = pstdd(i);% + sqrt(3)*pstdd(i);
      if para(i)>=a && para(i)<=b
          lnprior = lnprior + log(1/(b-a));
      else
          lnprior = lnprior + log(0);
      end
    end 
    
    prioinfo(i,1) = a;
    prioinfo(i,2) = b;
  end

  %
  %prioinfo;
  %
  % to run this part, introduce local variables x and y
  %graphset;
  %begwind;
  %window(3,3,0);
  %i = 1;
  %do until i > nprio;
  %   if prioinfo(i,3) == 1;
  %      x = seqa(0.01,0.01,98);
  %      y = exp(lpdfbeta(x,prioinfo(i,1),prioinfo(i,2)));
  %      xy(x,y);
  %   endif;
  %   nextwind;
  %   i=i+1;
  %endo;
  %endwind;

