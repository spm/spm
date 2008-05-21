function [dx] = spm_sde_dx(dfdx,dfdw,f,t)
% returns dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f + w for SDEs
% FORMAT [dx] = spm_sde_dx(dfdx,dfdw,f,t)
% dfdx   = df/dx - x: states
% dfdw   = df/dw - w: i.i.d. Weiner process 
% f      = dx/dt
% t      = integration time: (default t = 1);
%
% dx     = x(t) - x(0)
%--------------------------------------------------------------------------
% Integration of stochastic differential equations using local linearization. 
% This scheme accommodates nonlinearities in the state equation by using a 
% functional of f(x) = dx/dt.  This uses the equality
%
%             expm([0    0]*t) = expm(dfdx*t) - I)*inv(dfdx)*f
%                  [f dfdx]
%
% When t -> Inf this reduces to
%
%              dx(t) = -inv(dfdx)*f
%
% for the SDE:  dx = dfdx*x*dt + sqrt(2)*dfdw*dw
%
% where w is a standard Wiener process
%
% see also spm_dx
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_sde_dx.m 1703 2008-05-21 13:59:23Z karl $
 
% defaults
%--------------------------------------------------------------------------
if nargin < 3, t = 1; end
 
% compute stochastic terms {E = exp(dfdx*t), e = exp(dfdx*dt)}
%--------------------------------------------------------------------------
m     = length(dfdx);
N     = 256;
dt    = t/N;
eJdt  = expm(full(dfdx*dt));
eJt   = eJdt;
Q     = sparse(m,m);
R     = dfdw*dfdw'*2;
dQ    = eJt*R*eJt';
TOL   = trace(dQ)/64;
for i = 1:N
    
    % integrate and update exp(dfdx*t)
    %----------------------------------------------------------------------
    Q   = Q + dQ*dt;
    eJt = eJt*eJdt;
    dQ  = eJt*R*eJt';
    
    % convergence
    %----------------------------------------------------------------------
    if trace(dQ) < TOL, break, end
    
end
 
% scaled [Wiener] innovation
%--------------------------------------------------------------------------
w     = spm_sqrtm(Q)*randn(m,1);
 
% local linear solution plus stochastic term
%==========================================================================
dx    = spm_dx(dfdx,f,t) + w;
