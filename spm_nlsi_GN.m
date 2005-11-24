function [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using a Levenburg-Marquardt/EM algorithm
% FORMAT [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
% M.f  - f(x,u,P)
% M.g  - g(x,u,P)
%   x  - state  variables
%   u  - inputs or causes
%   P  - free parameters
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
% M.IS - integration scheme (function name or handle f(P,M,U,varargin))
%
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Ce - error precision constraints
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation  E{P|y}
% Cp  - (p x p)         conditional covariance   Cov{P|y}
% S   - (v x v)         [Re]ML estimate of       Cov{e}
%
% log evidence
%--------------------------------------------------------------------------
% F       - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% dynamic MIMO input-state-output model under Gaussian assumptions
%
%              dx/dt = f(x,u,P)
%              y     = g(x,u,P) + e
%
% A static nonlinear observation model with fixed input or causes u
% obtains when x = []. i.e.
%
%              y     = g(u,P) + e
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expectation and covariance of P
% If the free-energy starts to increase, a Levenburg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
%
%--------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_nlsi_GN.m 310 2005-11-24 16:27:02Z karl $

% figure
%--------------------------------------------------------------------------
Fsi = spm_figure;
set(Fsi,'name','System identification')
figure(Fsi)
 
 
% integration scheme
%--------------------------------------------------------------------------
try
    f = fcnchk(M.IS);
catch
    f = 'spm_int';
end
 
% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    M.x = sparse(M.n,1);
end
 
% data y
%--------------------------------------------------------------------------
y       = Y.y;
[ns nr] = size(y);          % number of samples and responses
 
% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Ce;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of precision components
ne    = length(Q{1});       % number of error terms
h     = zeros(nh,1);         % initialise hyperparameters
 
% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;
 
% confounds (if specified)
%--------------------------------------------------------------------------
try
    Ju = kron(speye(nr,nr),Y.X0);
catch
    Ju = sparse(ns*nr,0);
end
 
% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,1e-8);
nu    = size(Ju,2);                   % number of parameters (confounds)
np    = size(V, 2);                   % number of parameters (effective)
ip    = [1:np];
iu    = [1:nu] + np;
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
uC    = speye(nu)/1e-8;
ipC   = inv(spm_cat(diag({pC,uC})));
 
% initialise
%--------------------------------------------------------------------------
p     = sparse(np + nu,1);
Ep    = pE;
Cp    = pC;
 
C.F   = -Inf;
dv    = 1/128;
t     = 64;
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
hP    = eye(nh,nh)/16;
 
 
% EM
%==========================================================================
for k = 1:64
 
 
    % M-Step: ReML estimator of variance components:  h = max{F(p,h)}
    %======================================================================
 
    % prediction and errors
    %----------------------------------------------------------------------
    g     = feval(f,Ep,M,U,ns);
    e     = y(:) - g(:) - Ju*p(iu);
 
    % compute partial derivatives [Jp] df(p)/dp*V {p = parameters}
    %----------------------------------------------------------------------
    for i = ip
        dV      = dv*sqrt(Cp(i,i));
        dg      = feval(f,Ep + V(:,i)*dV,M,U,ns) - g;
        Jp(:,i) = dg(:)/dV;
    end
    J     = -[Jp Ju];
 
    % iterate a Fisher scoring scheme to find h = max{F(p,h)}
    %----------------------------------------------------------------------
    for m = 1:16
 
        % precision and conditional covariance
        %------------------------------------------------------------------
        iS    = speye(ne,ne)*1e-8;
        for i = 1:nh
            iS = iS + Q{i}*exp(h(i));
        end
        S     = inv(iS);
        Cp    = inv(J'*iS*J + ipC);
 
        % precision operators for M-Step
        %------------------------------------------------------------------
        for i = 1:nh
            P{i}  = Q{i}*exp(h(i));
            PS{i} = P{i}*S;
        end
 
        % derivatives: dLdh = dL/dh,...
        %------------------------------------------------------------------
        for i = 1:nh
            dFdh(i,1)      =  trace(PS{i})/2 - e'*P{i}*e/2 ...
                             -sum(sum(Cp.*(J'*P{i}*J)))/2;
            for j = i:nh
                dFdhh(i,j) = -sum(sum(PS{i}.*PS{j}))/2;
                dFdhh(j,i) =  dFdhh(i,j);
            end
        end
        
        % add hyperpriors
        %------------------------------------------------------------------
        dFdh  = dFdh  - hP*h;
        dFdhh = dFdhh - hP;
    
        % update ReML estimate
        %------------------------------------------------------------------
        Ch    = pinv(-dFdhh);
        dh    = Ch*dFdh;
        h     = h  + dh;
 
        % prevent overflow
        %------------------------------------------------------------------
        h     = max(h,-16);
        
        % convergence
        %------------------------------------------------------------------
        dF    = dFdh'*dh;
        if dF < 1e-2, break, end
 
    end
 
    % E-Step with Levenburg-Marquardt regularisation
    %======================================================================
 
    % objective function: F(p) (= log evidence - divergence)
    %----------------------------------------------------------------------
    F = - e'*iS*e/2 ...
        - p'*ipC*p/2 ...
        - ne*log(8*atan(1))/2 ...
        + spm_logdet(iS )/2 ...
        + spm_logdet(ipC)/2 ...
        + spm_logdet(Cp )/2 ...
        + spm_logdet(Ch )/2 ...
        - nh/2 - np/2;
 
 
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F
 
        % accept current estimates
        %------------------------------------------------------------------
        C.p   = p;
        C.h   = h;
        C.F   = F;
 
        % E-Step: Conditional update of gradients and curvature
        %------------------------------------------------------------------
        dFdp  = -J'*iS*e - ipC*p;
        dFdpp = -J'*iS*J - ipC;
 
        % decrease regularization
        %------------------------------------------------------------------
        t     = t*2;
        str   = 'E-Step(-)';
 
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        p     = C.p;
        h     = C.h;
 
        % and increase regularization
        %------------------------------------------------------------------
        t     = min(t/2,1);
        str   = 'E-Step(+)';
 
    end
 
    % E-Step: update
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{t});
    p     = p  + dp;
    Ep    = pE + V*p(ip);
    
    % graphics
    %----------------------------------------------------------------------
    figure(Fsi)

    % subplot parameters
    %----------------------------------------------------------------------
    subplot(2,1,1)
    plot([1:ns]*Y.dt,g),                        hold on
    plot([1:ns]*Y.dt,g + reshape(e,ns,nr),':'), hold off
    xlabel('time')
    title(sprintf('%s: %i','E-Step',k))
    grid on
    drawnow

    % subplot parameters
    %----------------------------------------------------------------------
    subplot(2,1,2)
    bar(full(V*p(ip)))
    xlabel('parameter')
    title('conditional [minus prior] expectation')
    grid on
    drawnow

    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    fprintf('%-6s: %i %6s %e %6s %e\n',str,k,'F:',C.F,'dp:',full(dF))
    if k > 2 && dF < 1e-2, break, end
    
end
 
% outputs
%--------------------------------------------------------------------------
Ep     = pE + V*p(ip);
Cp     = V*Cp(ip,ip)*V';
F      = C.F;



