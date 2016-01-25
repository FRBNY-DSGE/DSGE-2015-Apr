function [r,varargout] = distsmth_k93(y,pred,vpred,T,R,Q,Z,b,peachcount,psize,nant,antlags)

% DISTSMTH_K93.M

% This is a Kalman Smoothing program based on S.J. Koopman's "Disturbance
% Smoother for State Space Models" (Biometrika, 1993), as specified in
% Durbin and Koopman's "A Simple and Efficient Simulation Smoother for
% State Space Time Series Analysis" (Biometrika, 2002). The algorithm has been
% simplified for the case in which there is no measurement error, and the
% model matrices do not vary with time.

% This disturbance smoother is intended for use with the state smoother
% kalsmth_93.m from the same papers (Koopman 1993, Durbin and Koopman
% 2002). It produces a matrix of vectors, r, that is used for state 
% smoothing, and an optional matrix, eta_hat, containing the smoothed
% shocks. It has been adjusted to account for the possibility of missing
% values in the data, and to accommodate the zero bound model, which
% requires that no anticipated shocks occur before the zero bound window,
% which is achieved by setting the entries in the Q matrix corresponding to
% the anticipated shocks to zero in those periods.

% Nz will stand for the number of states, Ny for the number of observables,
% Ne for the number of shocks, and Nt for the number of periods of data.

% The state space is assumed to take the form:
% y(t) = Z*alpha(t) + b
% alpha(t+1) = T*alpha(t) + R*eta(t+1)

% INPUTS:

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
%       (i.e. the number of periods for which smoothed states are not
%       required).

% OUTPUTS:

% r, the (Nz x Nt) matrix used for state smoothing.
% eta_hat, the optional (Ne x Nt) matrix of smoothed shocks.

% Dan Greenwald, 7/7/2010.

if nargin < 10, error('At least ten inputs required'); end
if nargin == 11, error('Zero or two optional inputs required'); end
if nargin > 12, error('At most twelve inputs allowed'); end

Nt = size(y,2);
Nz = length(T);

r = zeros(Nz,Nt); % holds r_T-1,...r_0
r_t = zeros(Nz,1);

if nargout > 1
    Ne = size(R,2);
    eta_hat = zeros(Ne,Nt);
end

for t = Nt:-1:1
    
    y_t = y(:,t);
    
    % This section deals with the possibility of missing values in the y_t
    % vector (especially relevant for smoothing over peachdata).
    notnan = ~isnan(y_t);
    y_t = y_t(notnan);
    Z_t = Z(notnan,:);
    b_t = b(notnan);
    
    a = pred(:,t);
    P = vpred(:,:,t);
    
    F = Z_t*P*Z_t';
    v = y_t - Z_t*a - b_t;
    K = T*P*Z_t'/F;
    L = T - K*Z_t;
    
    r_t = Z_t'/F*v + L'*r_t;
    r(:,t) = r_t;
    
    if nargout > 1
        
        % This section relates to the zero bound framework, in which no
        % anticipated shocks are supposed to occur before the model switch.
        % In these periods, this is accomplished by setting the relevant
        % rows and columns of the Q matrix to zero. In other periods, or in
        % specifications with zero bound off (and hence with nant = 0), the
        % normal Q matrix can be used.
        
%         if exist('nant','var') && exist('antlags','var')
        if ~isempty(nant) && ~isempty(antlags)
            % The first part of the conditional below pertains to the periods in which zerobound is off. 
            % To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y. 
            % JC 11/30/10
            if nant > 0 && t < Nt-antlags-(peachcount*psize) 
                Q_t = zeros(Ne,Ne);
                Q_t(1:Ne-nant,1:Ne-nant) = Q(1:Ne-nant,1:Ne-nant);
                eta_hat(:,t) = Q_t*R'*r_t;
            else
                eta_hat(:,t) = Q*R'*r_t;
            end
        else
            eta_hat(:,t) = Q*R'*r_t; 
        end
    end
    
end

if nargout > 1, varargout(1) = {eta_hat}; end