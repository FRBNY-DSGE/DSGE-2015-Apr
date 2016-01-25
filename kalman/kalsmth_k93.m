function [alpha_hat,varargout] = kalsmth_k93(A0,P0,y,pred,vpred,T,R,Q,Z,b,nant,antlags,peachcount,psize,Ny0)

% KALSMTH_K93.M

% This is a Kalman Smoothing program based on S.J. Koopman's "Disturbance
% Smoother for State Space Models" (Biometrika, 1993), as specified in
% Durbin and Koopman's "A Simple and Efficient Simulation Smoother for
% State Space Time Series Analysis" (Biometrika, 2002). The algorithm has been
% simplified for the case in which there is no measurement error, and the
% model matrices do not vary with time.

% Unlike other Kalman Smoothing programs, there is no need to invert
% singular matrices using the Moore-Penrose pseudoinverse (pinv), which
% should lead to efficiency gains and fewer inversion problems. Also, the
% states vector and the corresponding matrices do not need to be augmented 
% to include the shock innovations. Instead they are saved automatically 
% in the eta_hat matrix.

% Nz will stand for the number of states, Ny for the number of observables,
% Ne for the number of shocks, and Nt for the number of periods of data.

% The state space is assumed to take the form:
% y(t) = Z*alpha(t) + b
% alpha(t+1) = T*alpha(t) + R*eta(t+1)

% INPUTS:

% A0, the (Nz x 1) initial (time 0) states vector.
% P0, the (Nz x Nz) initial (time 0) state covariance matrix.
% y, the (Ny x Nt) matrix of observable data.
% pred, the (Nz x Nt) matrix of one-step-ahead predicted states (from the Kalman Filter).
% vpred, the (Nz x Nz x Nt) matrix of one-step-ahead predicted covariance matrices.
% T, the (Nz x Nz) transition matrix.
% R, the (Nz x Ne) matrix translating shocks to states.
% Q, the (Ne x Ne) covariance matrix for the shocks.
% Z, the (Ny x Nz) measurement matrix.
% b, the (Ny x 1) constant vector in the measurement equation.

% nant, an optional scalar for the zero bound specification indicating the 
%       number of periods ahead the interest rate is fixed.
% antlags, an optional scalar for the zero bound specification indicating
%       the number of periods for which interest rate expectations have
%       been fixed
% Ny0, an optional scalar indicating the number of periods of presample
%       (i.e. the number of periods for which smoothed states are not required).

% OUTPUTS:

% alpha_hat, the (Nz x Nt) matrix of smoothed states.
% eta_hat, the optional (Ne x Nt) matrix of smoothed shocks.

% If Ny0 is nonzero, the alpha_hat and eta_hat matrices will be shorter by
% that number of columns (taken from the beginning).

% Dan Greenwald, 7/7/2010.

% if nargin < 12, error('At least twelve inputs required'); end
% if nargin == 13, error('Zero, two, or three optional inputs required'); end
% if nargin > 15, error('At most 15 inputs allowed'); end

Ne = size(R,2);
Nt = size(y,2);
% Ny = size(y,1);
Nz = length(T);

alpha_hat = zeros(Nz,Nt);

if nargout > 1
    [r,eta_hat] = distsmth_k93(y,pred,vpred,T,R,Q,Z,b,peachcount,psize,nant,antlags);
else
    r = distsmth_k93(y,pred,vpred,T,R,Q,Z,b,peachcount,psize);
end

ah_t = A0 + P0*r(:,1);
alpha_hat(:,1) = ah_t;

for t = 2:Nt

    % This section relates to the zero bound framework, in which no
    % anticipated shocks are supposed to occur before the model switch.
    % In these periods, this is accomplished by setting the relevant
    % rows and columns of the Q matrix to zero. In other periods, or in
    % specifications with zero bound off (and hence with nant = 0), the
    % normal Q matrix can be used.
    
    if ~isempty(antlags) && ~isempty(nant);
        % The first part of the conditional below pertains to the periods in which zerobound is off. 
        % To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y. 
        % JC 11/30/10
        if nant > 0 && t < Nt-antlags-(peachcount*psize) 
            Q_t = zeros(Ne,Ne);
            Q_t(1:Ne-nant,1:Ne-nant) = Q(1:Ne-nant,1:Ne-nant);
            ah_t = T*ah_t + R*Q_t*R'*r(:,t);
        else
            ah_t = T*ah_t + R*Q*R'*r(:,t);
        end
    else
        ah_t = T*ah_t + R*Q*R'*r(:,t);
    end
    
    alpha_hat(:,t) = ah_t;
end

if exist('Ny0','var')
   alpha_hat = alpha_hat(:,Ny0+1:end);
   if nargout > 1, eta_hat = eta_hat(:,Ny0+1:end); end
end

if nargout > 1, varargout(1) = {eta_hat}; end