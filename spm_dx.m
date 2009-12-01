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
% $Id: spm_dx.m 3605 2009-12-01 13:29:43Z karl $

% defaults
%--------------------------------------------------------------------------

if nargin < 3, t = Inf; end

% t is a regulariser
%--------------------------------------------------------------------------
sw = warning('off','MATLAB:log:logOfZero');
if iscell(t)
    t  = exp(t{:} - log(diag(-dfdx)));
end
warning(sw);

% use a [pseudo]inverse if all t > TOL
%==========================================================================
if min(t) > exp(16)

    dx = -spm_pinv(dfdx)*spm_vec(f);
    dx =  spm_unvec(dx,f);

else
    
    % ensure t is a scalar or matrix
    %----------------------------------------------------------------------
    if isvector(t), t = diag(t); end

    % augment Jacobian and take matrix exponential
    %======================================================================
    q  = sparse(1,1,1,length(dfdx) + 1,1);
    dx = spm_expm((spm_cat({0 []; t*spm_vec(f) t*dfdx})),q);
    dx = spm_unvec(dx(2:end),f);
    
end
dx     = real(dx);

