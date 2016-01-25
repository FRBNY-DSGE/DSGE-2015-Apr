function [hessian,stoph] = hessizero(fcn,x,Verbose,varargin)
% computes hessian of function fcn (string) evaluated at x (vector)
% varargin are the other inputs of fcn
% if Verbose, display error messages , results, etc.
% 11/12/01 translated by Marco DelNegro in matlab from Frank Schorfheide's program in gauss

%% index of free parameters
para_free = 1-x(:,2);
fpara_free = find(para_free);
nfree = length(fpara_free);

%% actual max
x = x(:,1);

npara = length(x);
ndx = 6;
dx =  exp(-(6:2:(6+(ndx-1)*2))');
hessian = zeros( npara, npara );
gradx = zeros(ndx,1);
grady = zeros(ndx,1);
gradxy = zeros(ndx,1);
hessdiag = zeros(ndx,1);
dxscale = ones(npara,1);



% Compute Diagonal elements first
for seli = fpara_free'

	if Verbose; fprintf(1,'\n Hessian Element: (%2.2g %2.2g)',[seli,seli]); end;
	for i=1:ndx;
		paradx = x;
		parady = x;
		paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
		parady(seli) = parady(seli) - dx(i)*dxscale(seli);
		paradxdy = paradx;
		paradxdy(seli) = paradxdy(seli) - dx(i)*dxscale(seli);
fx  = eval([fcn '(x,varargin{:})']);
		fdx = eval([fcn '(paradx,varargin{:})']);
		fdy = eval([fcn '(parady,varargin{:})']);
        fdxdy = eval([fcn '(paradxdy,varargin{:})']);
		gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
		grady(i) = ( fx - fdy )/ (dx(i)*dxscale(seli));
		gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(seli))^2 + (dx(i)*dxscale(seli))^2 );
		hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
		hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(seli));
    end
  	if Verbose == 2; fprintf(1,'\n Values: %2.6f',-hessdiag);         pause;     end;
	hessian(seli,seli) = -0.5*(hessdiag(3)+hessdiag(4));
    if hessian(seli,seli) <0
        error('negative diagonal in hessian');
    end
	if Verbose; fprintf(1,'\n Value Used: %2.6f ',hessian(seli,seli)); end;
end

% Now compute off-diagonal elements
% Make sure that correlations are between -1 and 1
% errorij contains the index of elements that are invalid
errorij = [ ];

for II = 1:(nfree-1);
   seli = fpara_free(II);	
   for JJ = II+1:nfree;
	selj = fpara_free(JJ);
    	if Verbose; fprintf(1,'\n Hessian Element: (%2.2g %2.2g)',[seli,selj]); end;
	    for i=1:ndx;
			paradx = x;
			parady = x;
			paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
			parady(selj) = parady(selj) - dx(i)*dxscale(selj);
			paradxdy = paradx;
			paradxdy(selj) = paradxdy(selj) - dx(i)*dxscale(selj);
    		fx  = eval([fcn '(x,varargin{:})']);
    		fdx = eval([fcn '(paradx,varargin{:})']);
    		fdy = eval([fcn '(parady,varargin{:})']);
            fdxdy = eval([fcn '(paradxdy,varargin{:})']);
			gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
			grady(i) = ( fx - fdy )/ (dx(i)*dxscale(selj));
			gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(selj))^2 + (dx(i)*dxscale(seli))^2 );
			hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
			hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(selj));
        end
    	if Verbose == 2; fprintf(1,'\n Values: %2.6f',-hessdiag); pause; end;
		%"Values"; 		;
		
		hessian(seli,selj) = -0.5*(hessdiag(3)+hessdiag(4));
		
		if ( hessian(seli,selj) == 0 ) | (hessian(selj,selj) == 0)
    		corrij = 0;
		else
		    corrij = hessian(seli,selj)/sqrt(hessian(seli,seli)*hessian(selj,selj));
        end
		
		if (corrij < -1) | (corrij > 1)
    		hessian(seli,selj)=0;
		    errorij = [errorij;[seli,selj,corrij]];
        end
		hessian(selj,seli) = hessian(seli,selj);
		
    	if Verbose 
            fprintf(1,'\n Value Used: %2.6f ',hessian(seli,seli)); 
            fprintf(1,'\n Correlation: %2.6f ',corrij); 
            fprintf(1,'\n Number of errors: %2.2g ',size(errorij,1)); 
        end
    end
end

	    stoph = 0;
if ~isempty(errorij)
    fprintf(1,'\n Errors: %2.6f %2.6f %2.6f  ',errorij');
    stoph = 1;
end


   
