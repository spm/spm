function [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using a Levenburg-Marquardt/EM algorithm
% FORMAT [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
% M.IS - function name or handle f(P,M,U)
%        This function specifies the nonlinear model: 
%        y = Y.y = IS(P,M,U) + X0*P0 + e
%        were e ~ N(0,C).  For dunamic systems this would be an intgration
%        scheme (e.g. spm_int). spm_int expects the following:
%
%     M.f  - f(x,u,P)
%     M.g  - g(x,u,P)
%       x  - state variables
%       u  - inputs or causes
%       P  - free parameters
%
% M.FS - function name or handle f(P,M,U,varargin)
%        This [optional] function perfoms feature selection assuming the
%        generlised model y = FS(Y,E) = FS(IS(P,M,U)) + e
%        where E are the parameters of the selection function
%
% M.E  - starting estimates for selection parameters
% M.P  - starting estimtes for model parameters [optional]
%
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
%
% M.hE - prior expectation  - E{h}
% M.hC - prior covariance   - Cov{h}
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Q  - error precision components
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
% nonliner model speficied by IS(P,M,U) under Gaussian assumptions. Usually,
% IS would be an integrator of a dynamic MIMO input-state-output model 
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
% An option feature slection can be specified with paramers M.E.
%
%--------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_nlsi_GN.m 681 2006-11-16 20:21:15Z karl $

% figure (unless disabled)
%--------------------------------------------------------------------------
try
    M.nograph;
catch
    if strcmp(get(gcf,'name'),'System identification')
        Fsi = gcf;
    else
        Fsi = spm_figure;
        set(Fsi,'name','System identification')
        figure(Fsi)
    end
end

% prediction scheme (usually an integrator)
%--------------------------------------------------------------------------
try
    IS = fcnchk(M.IS,'P','M','U');
catch
    IS = fcnchk('spm_int','P','M','U');
end

% feature slection
%--------------------------------------------------------------------------
try
    FS  = fcnchk(M.FS,'y','E');
    try
        E = M.E;
    catch
        E = [];
    end
catch
    FS  = inline('y','y','E');
    E   = [];
end

% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    if ~isfield(M,'n'), M.n = 0;    end
    M.x = sparse(M.n,1);
end

% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end

% initial parameters
%--------------------------------------------------------------------------
try
    spm_vec(M.P) - spm_vec(M.pE);
catch
    M.P = M.pE;
end

% time-step
%--------------------------------------------------------------------------
try
    Y.dt;
catch
    Y.dt = 1;
end

 
% data y
%--------------------------------------------------------------------------
y       = FS(Y.y,E);
[ns nr] = size(y);          % number of samples and responses
M.ns    = ns;               % store in M.ns for integrator

% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of precision components
ne    = length(Q{1});       % number of error terms
h     = zeros(nh,1);        % initialise hyperparameters
 
% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;
 
% confounds (if specified)
%--------------------------------------------------------------------------
try
    dgdu = kron(speye(nr,nr),Y.X0);
catch
    dgdu = sparse(ns*nr,0);
end

% hyperpriors - expectation
%--------------------------------------------------------------------------
try
    hE = M.hE;
catch
    hE = sparse(nh,1) - 4;
end

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    hP = inv(M.hC);
catch
    hP = eye(nh,nh)/4;
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,1e-8);
nu    = size(dgdu,2);                 % number of parameters (confounds)
np    = size(V,2);                    % number of parameters (effective)
ip    = [1:np];
iu    = [1:nu] + np;
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
uC    = speye(nu)/1e-8;
ipC   = inv(spm_cat(diag({pC,uC})));
 
% initialise conditional density
%--------------------------------------------------------------------------
p     = [V'*(spm_vec(M.P) - spm_vec(M.pE)); sparse(nu,1)];
Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);;
Cp    = pC;
 
% EM
%==========================================================================
C.F   = -Inf;
dv    = 1/128;
t     = 64;
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
for k = 1:128
 
 
    % M-Step: ReML estimator of variance components:  h = max{F(p,h)}
    %======================================================================
 
    % prediction (of features)
    %----------------------------------------------------------------------
    g     = feval(IS,Ep,M,U);
    g     = feval(FS,g,E);
    
    % compute partial derivatives; dgdp
    %----------------------------------------------------------------------
    for i = ip
        dV        = dv*sqrt(Cp(i,i));
        pi        = spm_unvec(spm_vec(Ep) + V(:,i)*dV,pE);
        gi        = feval(IS,pi,M,U);
        dg        = feval(FS,gi,E)  - g;
        dgdp(:,i) = spm_vec(dg)/dV;
    end

    % prediction error and gradients
    %----------------------------------------------------------------------
    e     = spm_vec(y) - spm_vec(g) - dgdu*p(iu);
    J     = -[dgdp dgdu];
 
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
        dFdh  = dFdh  - hP*(h - hE);
        dFdhh = dFdhh - hP;
    
        % update ReML estimate
        %------------------------------------------------------------------
        Ch    = inv(-dFdhh);
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
        + spm_logdet(Ch )/2;
 
 
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
        str   = 'EM-Step(-)';
 
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        p     = C.p;
        h     = C.h;
 
        % and increase regularization
        %------------------------------------------------------------------
        t     = min(t/2,1);
        str   = 'EM-Step(+)';
 
    end
 
    % E-Step: update
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{t});
    p     = p + dp;
    Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);

    % graphics
    %----------------------------------------------------------------------
    try
        figure(Fsi)

        % subplot prediction
        %----------------------------------------------------------------------
        subplot(2,1,1)
        plot([1:ns]*Y.dt,g),                      hold on
        plot([1:ns]*Y.dt,g + spm_unvec(e,g),':'), hold off
        xlabel('time')
        title(sprintf('%s: %i','E-Step',k))
        grid on

        % subplot parameters
        %----------------------------------------------------------------------
        subplot(2,1,2)
        bar(full(V*p(ip)))
        xlabel('parameter')
        title('conditional [minus prior] expectation')
        grid on
        drawnow
    end

    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    fprintf('%-6s: %i %6s %e %6s %e\n',str,k,'F:',C.F,'dF:',full(dF))
    if k > 2 && dF < 1e-2, break, end

end
 
% outputs
%--------------------------------------------------------------------------
Ep     = spm_unvec(spm_vec(pE) + V*p(ip),pE);;
Cp     = V*Cp(ip,ip)*V';
F      = C.F;



