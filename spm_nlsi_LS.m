function [Ep,qC,qh,F] = spm_nlsi_LS(M,U,Y)
% Bayesian inversion of a nonlinear model using (Laplacian) sampling
% FORMAT [Ep,Cp,Eh,F] = spm_nlsi_LS(M,U,Y)
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
% and covariance pC.
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
% $Id: spm_nlsi_LS.m 5219 2013-01-29 17:07:07Z spm $

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
h     = zeros(nh,1);                % initialise hyperparameters



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
        hE = hE*sparse(nh,1);
    end
catch
    hE = sparse(nh,1);
end

% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;
nu    = size(dfdu,2);                 % number of parameters (confounds)
np    = size(pC,2);                   % number of parameters (effective)

% second-order moments (in reduced space)
%--------------------------------------------------------------------------
ipC   = spm_inv(pC);


% initialize conditional density
%--------------------------------------------------------------------------
Eu    = spm_pinv(dfdu)*spm_vec(y);
Ep    = pE;

% precision and conditional covariance
%------------------------------------------------------------------
iS    = sparse(0);
for i = 1:nh
    iS = iS + Q{i}*(exp(-16) + exp(hE(i)));
end
S     = spm_inv(iS);
iS    = kron(speye(nq),iS);
qS    = spm_sqrtm(pC/32);
qE    = spm_vec(pE);
pE    = spm_vec(pE);
y     = spm_vec(y);
np    = length(qE);


% Sampling
%==========================================================================
Gmax  = -Inf;
for k = 1:64
    
    % time
    %----------------------------------------------------------------------
    tic;
    
    % Gibb's sampling
    %======================================================================
    for i = 1:128
        
        % prediction
        %------------------------------------------------------------------
        P(:,i) = qE + qS*randn(np,1);
        R(:,i) = spm_vec(feval(IS,spm_unvec(P(:,i),M.pE),M,U));

        % prediction error
        %------------------------------------------------------------------
        ey     = R(:,i) - y;
        ep     = P(:,i) - pE;
        
        % Gibb's energy
        %------------------------------------------------------------------
        qh     = real(ey')*iS*real(ey) + imag(ey)'*iS*imag(ey);
        G(i,1) = - ns*log(qh)/2 - ep'*ipC*ep/2;
        
        
        % conditional mode
        %----------------------------------------------------------------------
        [maxG,j] = max(G);
        if maxG  > Gmax
            qE   = P(:,j);
            f    = R(:,j);
            Gmax = maxG;
        end
        pE       = qE;
        
        disp(i)
    end
    
    
    % conditional dispersion
    %----------------------------------------------------------------------
    q     = exp((G - maxG));
    q     = q/sum(q);
    for i = 1:np
        for j = 1:np
            qC(i,j) = ((P(i,:) - qE(i)).*(P(j,:) - qE(j)))*q;
        end
    end
    qS    = spm_sqrtm(qC);
      
    
    % objective function:
    %======================================================================
    F     = Gmax + spm_logdet(ipC*qC)/2;
    F     = Gmax;
    
    
    % graphics
    %----------------------------------------------------------------------
    if exist('Fsi', 'var')
        spm_figure('Select', Fsi)
        
        % reshape prediction if necessary
        %------------------------------------------------------------------
        f  = reshape(f,ns,nr);
        d  = reshape(y,ns,nr);
        
        % subplot prediction
        %------------------------------------------------------------------
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
            plot(x,f,x,d,':')
            xlabel(xLab)
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
            
        else
            subplot(2,2,1)
            plot(x,real(f),x,real(d),':')
            xlabel(xLab)
            ylabel('real')
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
            
            subplot(2,2,2)
            plot(x,imag(f),x,imag(d),':')
            xlabel(xLab)
            ylabel('imaginary')
            title(sprintf('%s: %i','prediction and response: E-Step',k))
            grid on
        end
        
        % subplot Gibb's smapling
        %------------------------------------------------------------------
        subplot(2,2,3)
        plot(G)
        xlabel('smaple')
        title('Gibbs energy')
    
        % subplot parameters
        %------------------------------------------------------------------
        subplot(2,2,4)
        bar(full(qE - spm_vec(M.pE)))
        xlabel('parameter')
        title('conditional expectation')
        grid on
        drawnow
        
    end
    
    % convergence
    %----------------------------------------------------------------------
    try, dF = F - Fk; catch, dF = 0; end
    Fk      = F;
    fprintf('%-6s: %i %6s %-6.3e %6s %.3e (%.2f)\n','LS',k,'F:',full(F),'dF:',full(dF),toc)
    if k > 4 && dF < 1e-4
        break
    end  
end
if exist('Fsi', 'var')
    spm_figure('Focus', Fsi)
end
