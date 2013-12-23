function [Q,J] = spm_dcm_delay(M,P,D)
% returns the delay operator for flow and Jacobians of dynamical systems
% FORMAT [Q,J] = spm_dcm_delay(M,P,D)
% FORMAT [Q,J] = spm_dcm_delay(M,P)
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
% If the delay martix is not specifed it is computed from its parameters in
% P.D
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_delay.m 5817 2013-12-23 19:01:36Z karl $


% evaluate delay matrix D from parameters
%==========================================================================
if nargin < 3
    
    % paramterised delays
    %----------------------------------------------------------------------
    if isfield(P,'D')
        
        % number of states per sources
        %------------------------------------------------------------------
        nx  = size(M.x,2);
        
        % get prior means (log-delays)
        %------------------------------------------------------------------
        if isfield(M,'pF')
            di = M.pF.D(1);                    % intrinsic delays (ms)
            de = M.pF.D(2);                    % extrinsic delays (ms)
        else
            di = 1;
            de = 16;
        end
        
        % delay matrix D
        %------------------------------------------------------------------
        De  = exp(P.D);
        Di  = diag(diag(De));
        De  = De - Di;
        De  = De*de/1000;
        Di  = Di*di/1000;
        De  = kron(ones(nx,nx),De);
        Di  = kron(ones(nx,nx) - speye(nx,nx),Di);
        D   = Di + De;
        
    else
        
        % no delays
        %------------------------------------------------------------------
        if nargout < 2
            Q = sparse(1);
            return
        else
            D = sparse(0);
        end
    end
end

 
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
