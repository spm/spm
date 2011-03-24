function [Q,J] = spm_dcm_delay(M,P,D)
% returns the delay operator for flow and Jacobians of dynamical systems
% FORMAT [Q,J] = spm_dcm_delay(M,P,D)
%
% M   - model specification structure
% Required fields:
%   M.f - dx/dt    = f(x,u,P,M)            {function string or m-file}
%   M.m - m inputs
%   M.n - n states
%   M.x - (n x 1) = x(0) = expansion point: defaults to x = 0;
%   M.u - (m x 1) = u    = expansion point: defaults to u = 0;
%
% P     - model parameters
% D     - delay matrix (among hidden states): default D = zeros(n,n)
%
% return the delay operator for Jacobians of dynamical systems where the
% states are
%
% f     - dx(t)/dt  = f(x(t))
% Q     - delay operator dx(t)/dt = f(x(t - d))
%                                 = Q(d)*f(x(t))
% J     - Jacobian  = df/dt = (where delayed Jacobian = Q*J)
%
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_delay.m 4261 2011-03-24 16:39:42Z karl $
 
% create inline functions
%--------------------------------------------------------------------------
try
    funx = fcnchk(M.f,'x','u','P','M');
catch
    M.f  = inline('sparse(0,1)','x','u','P','M');
    M.n  = 0;
    M.x  = sparse(0,0);
    funx = fcnchk(M.f,'x','u','P','M');
end
 
% expansion point
%--------------------------------------------------------------------------
try, x = spm_vec(M.x); catch,  x = sparse(M.n,1); end
try, u = spm_vec(M.u); catch,  u = sparse(M.m,1); end
 
% Jacobian and delay operator
%==========================================================================
 
% derivatives
%--------------------------------------------------------------------------
J     = spm_diff(funx,x,u,P,M,1);
 
% delay operator
%--------------------------------------------------------------------------
D     = -D;
QJ    = (speye(length(J)) - D.*J)\J;
Q     = J;
for n = 1:128
    
    % n-th order Taylor term
    %----------------------------------------------------------------------
    dQ = ((D.^n).*J)*(QJ^n)/factorial(n);
    Q  = Q + dQ;
    
    % break if convergence
    %----------------------------------------------------------------------
    if norm(dQ,'inf') < 1e-6;     break, end

end
Q      = Q*spm_inv(J);
