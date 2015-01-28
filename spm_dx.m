function [dx] = spm_dx(dfdx,f,t)
% returns dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f
% FORMAT [dx] = spm_dx(dfdx,f,[t])
% dfdx   = df/dx
% f      = dx/dt
% t      = integration time: (default t = Inf);
%          if t is a cell (i.e., {t}) then t is set to:
%          exp(t - log(diag(-dfdx))
%
% dx     = x(t) - x(0)
%--------------------------------------------------------------------------
% Integration of a dynamic system using local linearization.  This scheme
% accommodates nonlinearities in the state equation by using a functional of
% f(x) = dx/dt.  This uses the equality
%
%             expm([0   0     ]) = (expm(t*dfdx) - I)*inv(dfdx)*f
%                  [t*f t*dfdx]
%
% When t -> Inf this reduces to
%
%              dx(t) = -inv(dfdx)*f
%
% These are the solutions to the gradient ascent ODE
%
%            dx/dt   = k*f = k*dfdx*x =>
%
%            dx(t)   = expm(t*k*dfdx)*x(0)
%                    = expm(t*k*dfdx)*inv(dfdx)*f(0) -
%                      expm(0*k*dfdx)*inv(dfdx)*f(0)
%
% When f = dF/dx (and dfdx = dF/dxdx), dx represents the update from a
% Gauss-Newton ascent on F.  This can be regularised by specifying {t}
% A heavy regularization corresponds to t = -4 and a light
% regularization would be t = 4. This version of spm_dx uses an augmented
% system and the Pade approximation to compute requisite matrix
% exponentials
%
% references:
%
% Friston K, Mattout J, Trujillo-Barreto N, Ashburner J, Penny W. (2007).
% Variational free energy and the Laplace approximation. NeuroImage.
% 34(1):220-34
%
% Ozaki T (1992) A bridge between nonlinear time-series models and
% nonlinear stochastic dynamical systems: A local linearization approach.
% Statistica Sin. 2:113-135.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dx.m 6321 2015-01-28 14:40:44Z karl $

% defaults
%--------------------------------------------------------------------------
nmax  = 512;                        % threshold for numerical approximation
if nargin < 3, t = Inf; end         % integration time
xf    = f; f = spm_vec(f);          % vectorise
n     = length(f);                  % dimensionality
U     = 1;                          % default (null) projector

% t is a regulariser
%--------------------------------------------------------------------------
sw  = warning('off','MATLAB:log:logOfZero');
SVD = 0;
if iscell(t)
    
    
    % work in natural gradients
    %----------------------------------------------------------------------
    if any(any(abs(dfdx - dfdx') > exp(-16))) && SVD
        U    = spm_svd(dfdx,0);
        dfdx = U'*dfdx*U;
        f    = U'*f;
    end
    
    % relative integration time
    %----------------------------------------------------------------------
    t      = t{:};
    if isscalar(t)
        t  = exp(t - spm_logdet(dfdx)/n);
    else
        t  = exp(t - log(diag(-dfdx)));
    end
    
end
warning(sw);

% use a [pseudo]inverse if all t > TOL
%==========================================================================
if min(t) > exp(16)

    dx = -spm_pinv(dfdx)*f;

elseif n > nmax && max(t) < 2 && isscalar(t)
    
    % numerical approximation (removing very dissipative modes)
    % fprintf('\nnumerical approximation: n = %i',length(f));
    %----------------------------------------------------------------------
    tol     = 1e-6;
    OPT.tol = tol*norm((dfdx),'inf');
      
    n     = 512;
    [V,D] = eigs(dfdx,1,'SR',OPT);
    N     = max(-real(D)*2,n);
    dt    = 1/N;
    
    s     = 0;
    x     = f - f;
    dx    = f*dt;
    df    = dfdx*dx;
    while s < t
        
        % while there are dissipative modes
        %------------------------------------------------------------------
        while N > n;
            
            
            % while principal mode dissipiates
            %--------------------------------------------------------------
            while abs(df'*V) > tol && s < t
                for i = 1:8
                    df = dfdx*dx;
                    f  = f + df;
                    dx = f*dt;
                    x  = x + dx;
                    s  = s + dt;
                end
            end
            
            % sliminate the principal mode 
            %--------------------------------------------------------------
            dfdx  = dfdx - V*(pinv(V)*dfdx);
            [V,D] = eigs(dfdx,1,'SR',OPT);
            N     = max(-real(D)*2,n);
            dt    = 1/N; 
           
        end
        
        % continue integration until s > t
        %------------------------------------------------------------------
        f  = f + dfdx*dx;
        dx = f*dt;
        x  = x + dx;
        s  = s + dt;
        
    end
    
    dx = U*x;
    
else
    
    % ensure t is a scalar or matrix
    %----------------------------------------------------------------------
    if isvector(t), t = diag(t); end

    % augment Jacobian and take matrix exponential
    %======================================================================
    J = full(spm_cat({0   [];
                      t*f t*dfdx}));
                  
    % solve using matrix expectation
    %----------------------------------------------------------------------
    if n  <= nmax
        dx = expm(J);
    else
        [V,D] = eig(J,'nobalance'); 
        dx    = V*diag(exp(diag(D)))/V;
    end
    
    % recover update
    %----------------------------------------------------------------------
    dx = dx(2:end,1);
    
end
dx     = spm_unvec(U*real(dx),xf);

