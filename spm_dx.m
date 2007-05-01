function [dx] = spm_dx(dfdx,f,t)
% returns dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f
% FORMAT [dx] = spm_dx(dfdx,f,[t])
% dfdx   = df/dx
% f      = dx/dt
% t      = integration time: (default t = Inf);
%          if t is a cell (i.e., {t}) then t is set to t{1}/norm(dfdx)
%
% dx     = x(t) - x(0)
%--------------------------------------------------------------------------
% Integration of a dynamic system using local linearization.  This scheme
% accommodates nonlinearities in the state equation by using a functional of
% f(x) = dx/dt.  This uses the equality
%
%             expm([0    0]*t) = expm(dfdx*t) - I)*inv(dfdx)*f
%                  [f dfdx]
%
% When t -> Inf this reduces to
%
%              dx(t) = -inv(dfdx)*f
%
% When f = dF/dx (and dfdx = dF/dxdx), dx represents the update from a
% Gauss-Newton ascent on F.  This can be regularised by specifying a finite
% t, A heavy regularization corresponds to t = 1/norm(dfdx) and a light
% regularization would be t = 32/norm(dfdx).  norm(dfdx) represents an upper
% bound on the rate of convergence (c.f., a Lyapunov exponent of the
% ascent)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_dx.m 253 2005-10-13 15:31:34Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3, t = Inf;              ; end
if iscell(t),  t = t{1}/normest(dfdx); end

% use a [pseudo]inverse if t > 1e8
%==========================================================================
if t > 1e8
    try
        dx = -inv(dfdx)*f;
    catch
        dx = -pinv(full(dfdx))*f;
    end
    return
end

% augment Jacobian and take matrix exponential
%==========================================================================
Jx    = spm_cat({0 []; spm_vec(f) dfdx});
dx    = spm_expm(Jx*t);
dx    = spm_unvec(dx(2:end,1),f);
