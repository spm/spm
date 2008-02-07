function [x,P] = spm_ekf(M,y)
% Extended Kalman Filtering for dynamic models
% FORMAT [x,P] = spm_ekf(M,y)
% M - model specification structure
% y - output or data (N x T)
%
% M(1).x                            % initial states
% M(1).f  = inline(f,'x','v','P')   % state equation
% M(1).g  = inline(g,'x','v','P')   % observer equation
% M(1).pE                           % parameters
% M(1).V                            % observation noise precision
%
% M(2).v                            % initial process f(noise)
% M(2).V                            % process f(noise) precision
%
% x - conditional expectation of states
% P - {1 x T} conditional covariance of states
%__________________________________________________________________________
% See notes at the end of this script for details and a demo.  This routine
% is based on:
%
% var der Merwe R, Doucet A, de Freitas N and Wan E (2000). The
% unscented particle filter.  Technical Report CUED/F-INFENG/TR 380
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ekf.m 1143 2008-02-07 19:33:33Z spm $

% check model specification
%--------------------------------------------------------------------------
M  = spm_DEM_M_set(M);
dt = M(1).E.dt;
if length(M) ~=2
    errordlg('spm_ekf requires a two-level model')
    return
end

% INITIALISATION:
% =========================================================================
dfdx  = spm_diff(M(1).f,M(1).x,M(2).v,M(1).pE,1);
dfdv  = spm_diff(M(1).f,M(1).x,M(2).v,M(1).pE,2);
dgdx  = spm_diff(M(1).g,M(1).x,M(2).v,M(1).pE,1);
Jx    = spm_expm(dfdx);
T     = length(y);              % number of time points

% covariances
%--------------------------------------------------------------------------
iR    = M(1).V;
for i = 1:length(M(1).Q)
   iR = iR + M(1).Q{i}*exp(M(1).hE(i));
end
iQ    = M(1).W;
for i = 1:length(M(1).R)
   iQ = iQ + M(1).R{i}*exp(M(1).gE(i));
end
R  = inv(iR);                            % EKF measurement noise variance.
Q  = inv(iQ);                            % EKF process noise variance.
x  = M(1).x;                             % EKF estimate of the mean of states
Q  = Q + Jx*dfdv*inv(M(2).V)*dfdv'*Jx';  % EKF process noise variance.
P  = {pinv(full(dgdx'*R*dgdx))};         % EKF conditional covariance of states

for t = 2:T

    % PREDICTION STEP:
    %----------------------------------------------------------------------
    f        = M(1).f(M(1).x,M(2).v,M(1).pE);
    dfdx     = spm_diff(M(1).f,M(1).x,M(2).v,M(1).pE,1);
    xPred    = M(1).x + spm_dx(dfdx,f,dt);
    Jx       = spm_expm(dfdx);
    PPred    = Q + Jx*P{t-1}*Jx';

    % CORRECTION STEP:
    %----------------------------------------------------------------------
    yPred    = M(1).g(xPred,M(2).v,M(1).pE);
    Jy       = spm_diff(M(1).g,xPred,M(2).v,M(1).pE,1);
    K        = PPred*Jy'*inv(R + Jy*PPred*Jy');
    M(1).x   = xPred + K*(y(:,t) - yPred);
    x(:,t)   = M(1).x;
    P{t}     = PPred - K*Jy*PPred;

    % report
    %----------------------------------------------------------------------
    fprintf('EKF: time-step = %i : %i\n',t,T);

end


return

% notes and demo:
%==========================================================================
% The code below generates a nonlinear, non-Gaussian problem (S) comprising
% a model S.M and data S.Y (c.f. van der Merwe et al 2000))
%
% The model is   f(x) = dxdt
%                     = 1 + sin(o.o4*pi*t) - log(2)*x + n
%                y    = g(x)
%                     = (x.^2)/5  : if t < 30
%                       -2 + x/2  : otherwise
% i.e. the output nonlinearity becomes linear after 30 time steps.  In this
% implementation time is modelled as an auxiliary state variable.  n is
% the process noise, which is modelled as a log-normal variate.  e is
% Gaussian observation noise.

% model specification
%--------------------------------------------------------------------------
f       = '[1; (1 + sin(P(2)*pi*x(1)) - P(1)*x(2) + exp(v))]';
g       = '(x(1) > 30)*(-2 + x(2)/2) + ~(x(1) > 30)*(x(2).^2)/5';
M(1).x  = [1; 1];                  % initial states
M(1).f  = inline(f,'x','v','P');   % state equation
M(1).g  = inline(g,'x','v','P');   % observer equation
M(1).pE = [log(2) 0.04];           % parameters
M(1).V  = 1e5;                     % observation noise precision

M(2).v  = 0;                       % initial process log(noise)
M(2).V  = 2.4;                     % process log(noise) precision

% generate data (output)
%--------------------------------------------------------------------------
T       = 60;                      % number of time points
S       = spm_DEM_generate(M,T);

% EKF
%--------------------------------------------------------------------------
ekf_x   = spm_ekf(M,S.Y);


% plot results
%--------------------------------------------------------------------------
x       = S.pU.x{1};
plot([1:T],x(2,:),[1:T],ekf_x(2,:))
legend({'true','EKF'})
