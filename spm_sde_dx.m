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
% Integration of stochastic differential euquation using local linearization. 
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm__sde_dx.m 253 2005-10-13 15:31:34Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3, t = 1; end

% compute stochastic terms {E = exp(dfdx*t), e = exp(dfdx*dt)}
%--------------------------------------------------------------------------
m     = length(dfdx);
N     = 256;
dt    = t/N;
e     = spm_expm(dfdx*dt);
E     = e;
Q     = sparse(m,m);
R     = dfdw*dfdw'*2;
TOL   = trace(E*R*E')/64;
for i = 1:N
    Q = Q + E*R*E'*dt;
    E = E*e;
    if trace(E*R*E') < TOL, break, end
end
w     = spm_sqrtm(Q)*randn(m,1);

% local linear solution and add stochastic term
%==========================================================================
dx    = spm_dx(dfdx,f,t) + w;




