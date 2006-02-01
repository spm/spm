function [qx,qP] = spm_pf(M,y)
% Particle Filtering for dynamic models
% FORMAT [x,P] = spm_pf(M,y)
% M - model specification structure
% y - output or data (N x T)
%
% M(1).x                            % initial states
% M(1).f  = inline(f,'x','v','P')   % state equation
% M(1).g  = inline(g,'x','v','P')   % observer equation
% M(1).pE                           % parameters
% M(1).V                            % observation noise precision
%
% M(2).v                            % initial process noise
% M(2).V                            % process noise precision
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_pf.m 417 2006-02-01 13:50:14Z karl $

% check model specification
%--------------------------------------------------------------------------
M  = spm_M_set(M);
dt = M(1).E.dt;
if length(M) ~=2
    errordlg('spm_pf requires a two-level model')
    return
end

% INITIALISATION:
%==========================================================================
T  = length(y);                    % number of time points
N  = 200;                          % number of particles.

P  = M(1).P;                       % parameters
R  = M(1).V;                       % Precision of measurement noise
Q  = sqrtm(inv(full(M(2).V)));     % sqrt covariance of process noise
v  = M(2).v;
x  = kron(ones(1,N),M(1).x);

for t = 1:T

    % PREDICTION STEP: with the transition prior as proposal
    %----------------------------------------------------------------------
    for i = 1:N
        v          = Q*randn(size(v));
        f          = M(1).f(x(:,i),v,P);
        dfdx       = spm_diff(M(1).f,x(:,i),v,P,1);
        xPred(:,i) = x(:,i) + spm_dx(dfdx,f,dt);
    end

    % EVALUATE IMPORTANCE WEIGHTS: and normalise
    %----------------------------------------------------------------------
    for i = 1:N
        yPred  = M(1).g(xPred(:,i),M(2).v,P);
        ePred  = yPred - y(:,t);
        w(i)   = ePred'*R*ePred;
    end
    w   = w - min(w);
    w   = exp(-w/2);
    w   = w/sum(w);

    % SELECTION STEP: multinomial resampling.
    %----------------------------------------------------------------------    
    x   = xPred(:,multinomial(1:N,w));

    % report and record moments
    %----------------------------------------------------------------------
    qx(:,t)  = mean(x,2);
    qP{t}    = cov(x');
    fprintf('PF: time-step = %i : %i\n',t,T);

end

function I = multinomial(inIndex,q);
%==========================================================================
% PURPOSE : Performs the resampling stage of the SIR
%           in order(number of samples) steps.
% INPUTS  : - inIndex = Input particle indices.
%           - q       = Normalised importance ratios.
% OUTPUTS : - I = Resampled indices.
% AUTHORS : Arnaud Doucet and Nando de Freitas

% MULTINOMIAL SAMPLING:
% generate S ordered random variables uniformly distributed in [0,1]
% high speed Niclas Bergman Procedure
%--------------------------------------------------------------------------
q        = q(:);
S        = length(q);  % S = Number of particles.
N_babies = zeros(1,S);
cumDist  = cumsum(q');

u = fliplr(cumprod(rand(1,S).^(1./(S:-1:1))));
j = 1;
for i = 1:S
    while (u(1,i) > cumDist(1,j))
        j = j + 1;
    end
    N_babies(1,j) = N_babies(1,j) + 1;
end;

% COPY RESAMPLED TRAJECTORIES:
%--------------------------------------------------------------------------
index = 1;
for i = 1:S
    if (N_babies(1,i)>0)
        for j=index:index+N_babies(1,i)-1
            I(j) = inIndex(i);
        end;
    end;
    index = index + N_babies(1,i);
end

return
%==========================================================================

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
S       = spm_DEM_generate(M,sparse(1,T + 1));

% Particle filtering
%--------------------------------------------------------------------------
pf_x    = spm_pf(M,S.Y);

% plot results
%--------------------------------------------------------------------------
x       = S.pU.x{1};
plot([1:T],x(2,:),[1:T],pf_x(2,:))
legend({'true','PF'})
