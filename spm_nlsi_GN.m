function [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y)
% Bayesian inversion of a nonlinear model using a Gauss-Newton/EM algorithm
% FORMAT [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
%
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
% M.hE - prior expectation  - E{h}   of log-precision parameters
% M.hC - prior covariance   - Cov{h} of log-precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.dt - sampling interval for outputs
% Y.X0 - Confounds or null space      (over size(y,1) bins or all vec(y))
% Y.Q  - q error precision components (over size(y,1) bins or all vec(y))
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation    E{P|y}
% Cp  - (p x p)         conditional covariance     Cov{P|y}
% Eh  - (q x 1)         conditional log-precisions E{h|y}
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
% An optional feature selection can be specified with parameters M.FS.
%
% For generic aspects of the scheme see:
% 
% Friston K, Mattout J, Trujillo-Barreto N, Ashburner J, Penny W. 
% Variational free energy and the Laplace approximation. 
% NeuroImage. 2007 Jan 1;34(1):220-34.
% 
% This scheme handels complex data along the lines originally described in:
% 
% Sehpard RJ, Lordan BP, and Grant EH. 
% Least squares analysis of complex data with applications to permittivity
% measurements.
% J. Phys. D. Appl. Phys 1970 3:1759-1764.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nlsi_GN.m 4261 2011-03-24 16:39:42Z karl $
 
% figure (unless disabled)
%--------------------------------------------------------------------------
try
    M.nograph;
catch 
    M.nograph = 0;
end
if ~M.nograph
    Fsi = spm_figure('GetWin','SI');
end

% check integrator
%--------------------------------------------------------------------------
try
    M.IS;
catch
    M.IS = 'spm_int';
end
 
% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
try
    
    % try FS(y,M)
    %----------------------------------------------------------------------
    try
        y  = feval(M.FS,Y.y,M);
        IS = inline([M.FS '(' M.IS '(P,M,U),M)'],'P','M','U');
        
    % try FS(y)
    %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        IS = inline([M.FS '(' M.IS '(P,M,U))'],'P','M','U');
    end
    
catch
    
    % otherwise FS(y) = y
    %----------------------------------------------------------------------
    y   = Y.y;
    IS  = inline([M.IS '(P,M,U)'],'P','M','U');
end
 
% size of data (usually samples x channels)
%--------------------------------------------------------------------------
if iscell(y)
    ns = size(y{1},1);
else
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
    fprintf('\nParameter initialisation successful\n')
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
nq    = nr*ns/nt;                   % for compact Kronecker form of M-step
h     = sparse(nh,1);                % initialise hyperparameters
 
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
    if length(hE) ~= nh
        hE = hE(1) + sparse(nh,1);
    end
catch
    hE = sparse(nh,1);
end
h      = hE;
 
% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = inv(M.hC);
    if length(ihC) ~= nh
        ihC = ihC*speye(nh,nh);
    end
catch
    ihC = speye(nh,nh);
end
 
% unpack covariance
%--------------------------------------------------------------------------
if isstruct(pC);
    pC = spm_diag(spm_vec(pC));
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,exp(-32));
nu    = size(dfdu,2);                 % number of parameters (confounds)
np    = size(V,2);                    % number of parameters (effective)
ip    = (1:np)';
iu    = (1:nu)' + np;
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
uC    = speye(nu)/1e-8;
ipC   = inv(spm_cat(spm_diag({pC,uC})));
 
% initialize conditional density
%--------------------------------------------------------------------------
Eu    = spm_pinv(dfdu)*spm_vec(y);
p     = [V'*(spm_vec(M.P) - spm_vec(M.pE)); Eu];
Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);

 
% EM
%==========================================================================
criterion = [0 0 0 0];

C.F   = -Inf;                                   % free energy
v     = -2;                                     % log ascent rate
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
for k = 1:128
    
    % time
    %----------------------------------------------------------------------  
    tStart = tic;
 
    % M-Step: ReML estimator of variance components:  h = max{F(p,h)}
    %======================================================================
 
    % prediction f, and gradients; dfdp
    %----------------------------------------------------------------------
    [dfdp f] = spm_diff(IS,Ep,M,U,1,{V});
    
       
    % prediction error and full gradients
    %----------------------------------------------------------------------
    e     =  spm_vec(y) - spm_vec(f) - dfdu*p(iu);
    dfdp  =  reshape(spm_vec(dfdp),ns*nr,np);
    J     = -[dfdp dfdu];
    
 
    % M-step; Fisher scoring scheme to find h = max{F(p,h)}
    %======================================================================
    for m = 1:8
 
        % check for stability
        %------------------------------------------------------------------
        if norm(J,'inf') > exp(32), break, end
        
        % precision and conditional covariance
        %------------------------------------------------------------------
        iS    = sparse(0);
        for i = 1:nh
            iS = iS + Q{i}*(exp(-32) + exp(h(i)));
        end
        S     = spm_inv(iS);
        iS    = kron(speye(nq),iS);
        Pp    = real(J)'*iS*real(J) + imag(J)'*iS*imag(J);
        Cp    = spm_inv(Pp + ipC);
        if any(isnan(Cp(:))) || rcond(full(Cp)) < exp(-32), break, end
 
        % precision operators for M-Step
        %------------------------------------------------------------------
        for i = 1:nh
            P{i}   = Q{i}*exp(h(i));
            PS{i}  = P{i}*S;
            P{i}   = kron(speye(nq),P{i});
            JPJ{i} = real(J)'*P{i}*real(J) + ...
                     imag(J)'*P{i}*imag(J);
        end
                    
        % derivatives: dLdh = dL/dh,...
        %------------------------------------------------------------------
        for i = 1:nh
            dFdh(i,1)      =   trace(PS{i})*nq/2 ...
                             - real(e)'*P{i}*real(e)/2 ...
                             - imag(e)'*P{i}*imag(e)/2 ...
                             - sum(sum(Cp.*JPJ{i}))/2;
            for j = i:nh
                dFdhh(i,j) = - sum(sum(PS{i}.*PS{j}))*nq/2;
                dFdhh(j,i) =   dFdhh(i,j);
            end
        end
        
        % add hyperpriors
        %------------------------------------------------------------------
        d     = h     - hE;
        dFdh  = dFdh  - ihC*d;
        dFdhh = dFdhh - ihC;
        Ch    = spm_inv(-dFdhh); 
        
        % update ReML estimate
        %------------------------------------------------------------------
        dh    = spm_dx(dFdhh,dFdh,{4});
        dh    = min(max(dh,-1),1);
        h     = h  + dh;
 
        % convergence
        %------------------------------------------------------------------
        dF    = dFdh'*dh;
        if dF < 1e-2, break, end
 
    end

    
    % E-Step with Levenberg-Marquardt regularization
    %======================================================================
 
    % objective function: F(p) (= log evidence - divergence)
    %----------------------------------------------------------------------
    F = - real(e)'*iS*real(e)/2 ...
        - imag(e)'*iS*imag(e)/2 ...
        - p'*ipC*p/2 ...
        - d'*ihC*d/2 ...
        - ns*nr*log(8*atan(1))/2 ...
        - spm_logdet(S)*nq/2 ...
        + spm_logdet(ipC*Cp)/2 ...
        + spm_logdet(ihC*Ch)/2;
 
    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    try
        F0; 
        fprintf(' actual: %.3e (%.2f sec)\n',full(F - C.F),toc(tStart))
    catch
        F0 = F;
    end
 
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F
 
        % accept current estimates
        %------------------------------------------------------------------
        C.p   = p;
        C.h   = h;
        C.F   = F;
        C.Cp  = Cp;
        
        % E-Step: Conditional update of gradients and curvature
        %------------------------------------------------------------------
        dFdp  = -real(J)'*iS*real(e) ...
                -imag(J)'*iS*imag(e) - ipC*p;
        dFdpp = -real(J)'*iS*real(J) ...
                -imag(J)'*iS*imag(J) - ipC;
        
        % decrease regularization
        %------------------------------------------------------------------
        v     = min(v + 1/2,4);
        str   = 'EM:(+)';
        
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        p     = C.p;
        h     = C.h;
        Cp    = C.Cp;
 
        % and increase regularization
        %------------------------------------------------------------------
        v     = min(v - 2,-4);
        str   = 'EM:(-)';
        
    end
 
    % E-Step: update
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{v});
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
        x    = (1:ns)*Y.dt;
        xLab = 'time (seconds)';
        try
            if length(M.Hz) == ns
                x    = Y.Hz;
                xLab = 'Frequency (Hz)';
            end
        end
        
        if isreal(f)
            
            subplot(2,1,1)
            plot(x,f), hold on
            plot(x,f + spm_unvec(e,f),':'), hold off
            xlabel(xLab)
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
            
        else
            
            subplot(2,2,1)
            plot(x,real(f)), hold on
            plot(x,real(f + spm_unvec(e,f)),':'), hold off
            xlabel(xLab)
            ylabel('real')
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
            
            subplot(2,2,2)
            plot(x,imag(f)), hold on
            plot(x,imag(f + spm_unvec(e,f)),':'), hold off
            xlabel(xLab)
            ylabel('imaginary')
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
            
        end
        
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
    fprintf('%-6s: %i %6s %-6.3e %6s %.3e ',str,k,'F:',full(C.F - F0),'dF predicted:',full(dF))
    
    criterion = [(dF < 1e-2) criterion(1:end - 1)];
    if all(criterion), fprintf(' convergence\n'), break, end
 
end
 
% outputs
%--------------------------------------------------------------------------
Ep     = spm_unvec(spm_vec(pE) + V*C.p(ip),pE);
Cp     = V*C.Cp(ip,ip)*V';
Eh     = C.h;
F      = C.F;

