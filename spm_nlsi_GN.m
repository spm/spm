function [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
% Bayesian inversion of a nonlinear model using a Gauss-Newton/EM algorithm
% FORMAT [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
% M.IS - function name f(P,M,U) - generative model
%        This function specifies the nonlinear model: 
%        y = Y.y = IS(P,M,U) + X0*P0 + e
%        were e ~ N(0,C).  For dynamic systems this would be an integration
%        scheme (e.g. spm_int). spm_int expects the following:
%
%     M.f  - f(x,u,P,M)
%     M.g  - g(x,u,P,M)
%       x  - state variables
%       u  - inputs or causes
%       P  - free parameters
%       M  - fixed functional forms and parameters in M
%
% M.FS - function name f(y,M)   - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(IS(P,M,U),M) + X0*P0 + e
%
% M.P  - starting estimates for model parameters [optional]
%
% M.pE - prior expectation  - E{P}   of model parameters
% M.pC - prior covariance   - Cov{P} of model parameters
%
% M.hE - prior expectation  - E{h}   of precision parameters
% M.hC - prior covariance   - Cov{h} of precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.dt - sampling interval for outputs
% Y.X0 - Confounds or null space    (over size(y,1) bins or all vec(y))
% Y.Q  - error precision components (over size(y,1) bins or all vec(y))
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
% F   - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model specified by IS(P,M,U) under Gaussian assumptions. 
% Usually, IS is an integrator of a dynamic MIMO input-state-output model 
%
%              dx/dt = f(x,u,P)
%              y     = g(x,u,P)  + X0*P0 + e
%
% A static nonlinear observation model with fixed input or causes u
% obtains when x = []. i.e.
%
%              y     = g([],u,P) + X0*P0e + e
%
% but static nonlinear models are specified more simply using
%
%              y     = IS(P,M,U) + X0*P0 + e
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expectation and covariance of P
% If the free-energy starts to increase, a Levenberg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
%
% An optional feature selection can be specified with parameters M.FS
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nlsi_GN.m 1068 2008-01-07 18:53:03Z karl $
 
% figure (unless disabled)
%--------------------------------------------------------------------------
try
    M.nograph;
catch 
    Fsi = spm_figure('GetWin','SI');
end
 
% check integrator
%--------------------------------------------------------------------------
try
    M.IS;
catch
    M.IS = 'spm_int_U';
end
 
% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
if isfield(M,'FS')
    
    % FS(y,M)
    %----------------------------------------------------------------------
    try
        y  = feval(M.FS,Y.y,M);
        IS = inline([M.FS '(' M.IS '(P,M,U),M)'],'P','M','U');
        
    % FS(y,M)
    %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        IS = inline([M.FS '(' M.IS '(P,M,U))'],'P','M','U');
 
    end
else
    
    % FS(y) = y
    %----------------------------------------------------------------------
    y   = Y.y;
    IS  = inline([M.IS '(P,M,U)'],'P','M','U');
end
 
% size of data (usually samples x channels)
%--------------------------------------------------------------------------
try
    ns = M.ns;
catch
    ns = size(y,1);
end
nr   = length(spm_vec(y))/ns;       % number of samples and responses
M.ns = ns;                          % store in M.ns for integrator
 
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
 
 
% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
    if isnumeric(Q), Q = {Q}; end
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);                  % number of precision components
nt    = length(Q{1});               % number of time bins
nq    = nr*ns/nt;                   % for compact kronecker form of M-step
h     = zeros(nh,1);                % initialise hyperparameters
 
% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;
 
% confounds (if specified)
%--------------------------------------------------------------------------
try
    nb   = size(Y.X0,1);            % number of bins
    nx   = nr*ns/nb;                % number of blocks
    dfdu = kron(speye(nx,nx),Y.X0);
catch
    dfdu = sparse(ns*nr,0);
end
 
% hyperpriors - expectation
%--------------------------------------------------------------------------
try
    hE = M.hE;
catch
    hE = sparse(nh,1);
end
 
% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = inv(M.hC);
catch
    ihC = eye(nh,nh);
end
 
% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,exp(-16));
nu    = size(dfdu,2);                 % number of parameters (confounds)
np    = size(V,2);                    % number of parameters (effective)
ip    = [1:np];
iu    = [1:nu] + np;
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
uC    = speye(nu)/1e-8;
ipC   = inv(spm_cat(diag({pC,uC})));
 
% initialize conditional density
%--------------------------------------------------------------------------
Eu    = inv(dfdu'*dfdu)*(dfdu'*spm_vec(y));
p     = [V'*(spm_vec(M.P) - spm_vec(M.pE)); Eu];
Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);;
Cp    = pC;
 
% EM
%==========================================================================
C.F   = -Inf;
t     = 256;
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
warning off
for k = 1:128
 
 
    % M-Step: ReML estimator of variance components:  h = max{F(p,h)}
    %======================================================================
 
    % prediction f, and gradients; dfdp
    %----------------------------------------------------------------------
    [dfdp f] = spm_diff(IS,Ep,M,U,1,{V});
       
    % prediction error and full gradients
    %----------------------------------------------------------------------
    e     = spm_vec(y) - spm_vec(f) - dfdu*p(iu);
    dfdp  = reshape(spm_vec(dfdp),ns*nr,np);
    J     = -[dfdp dfdu];
    
 
    % M-step; Fisher scoring scheme to find h = max{F(p,h)}
    %======================================================================
    for m = 1:16
 
        % precision and conditional covariance
        %------------------------------------------------------------------
        iS    = speye(nt,nt)*1e-8;
        for i = 1:nh
            iS = iS + Q{i}*exp(h(i));
        end
        S     = inv(iS);
        iS    = kron(speye(nq),iS);
        Cp    = inv(J'*iS*J + ipC);
        
        
 
        % precision operators for M-Step
        %------------------------------------------------------------------
        for i = 1:nh
            P{i}  = Q{i}*exp(h(i));
            PS{i} = P{i}*S;
            P{i}  = kron(speye(nq),P{i});
        end
            
            
        % derivatives: dLdh = dL/dh,...
        %------------------------------------------------------------------
        for i = 1:nh
            dFdh(i,1)      =  trace(PS{i})*nq/2 - e'*P{i}*e/2 ...
                             -sum(sum(Cp.*(J'*P{i}*J)))/2;
            for j = i:nh
                dFdhh(i,j) = -sum(sum(PS{i}.*PS{j}))*nq/2;
                dFdhh(j,i) =  dFdhh(i,j);
            end
        end
        
        % add hyperpriors
        %------------------------------------------------------------------
        d     = h - hE;
        dFdh  = dFdh  - ihC*d;
        dFdhh = dFdhh - ihC;
    
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
 
    % E-Step with Levenberg-Marquardt regularization
    %======================================================================
 
    % objective function: F(p) (= log evidence - divergence)
    %----------------------------------------------------------------------
    F = - e'*iS*e/2 ...
        - p'*ipC*p/2 ...
        - d'*ihC*d/2 ...
        - ns*nr*log(8*atan(1))/2 ...
        - spm_logdet(S)*nq/2 ...
        + spm_logdet(ipC*Cp)/2 ...
        + spm_logdet(ihC*Ch)/2;
 
 
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
        t     = min(t/2,128);
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
 
        % reshape prediction if necessary
        %------------------------------------------------------------------
        f  = reshape(spm_vec(f),ns,nr);
        
        % subplot prediction
        %------------------------------------------------------------------
        figure(Fsi)
        subplot(2,1,1)
        plot([1:ns]*Y.dt,f),                      hold on
        plot([1:ns]*Y.dt,f + spm_unvec(e,f),':'), hold off
        xlabel('time')
        title(sprintf('%s: %i','E-Step',k))
        grid on
 
        % subplot parameters
        %------------------------------------------------------------------
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
warning on
