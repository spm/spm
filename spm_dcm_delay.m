function [Q,J] = spm_dcm_delay(P,M)
% returns the delay operator for flow and Jacobians of dynamical systems
% FORMAT [Q,J] = spm_dcm_delay(P,M)
% P   - model parameters
% M   - model specification structure
% Required fields:
%   M.f - dx/dt    = f(x,u,P,M)            {function string or m-file}
%   M.m - m inputs
%   M.n - n states
%   M.x - (n x 1) = x(0) = expansion point: defaults to x = 0;
%   M.u - (m x 1) = u    = expansion point: defaults to u = 0;
%
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
% P.D (and M.pF.D if specified)
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_delay.m 5964 2014-04-20 09:48:58Z karl $


% evaluate delay matrix D from parameters
%==========================================================================

% paramterised delays
%--------------------------------------------------------------------------
if isfield(P,'D')
    
    % number of states per sources
    %----------------------------------------------------------------------
    nx  = size(M.x,2);
    
    % get prior means (log-delays)
    %----------------------------------------------------------------------
    if isfield(M,'pF')
        di = M.pF.D(1);                    % intrinsic delays (ms)
        de = M.pF.D(2);                    % extrinsic delays (ms)
    else
        di = 1;
        de = 8;
    end
    
    % delay matrix D
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    D = sparse(0);
    
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
J     = full(spm_diff(funx,x,u,P,M,1));

% delay operator:  estimated using a Robbins–Monro algorithm
%--------------------------------------------------------------------------
N     = 256;
D     = -D;
QJ    = (eye(length(J)) - D.*J)\J;
a     = 1/2;
TOL   = norm(QJ,'inf')*1e-6;
dn    = 1;
Dn    = cell(N,1);
for n = 1:N
    dn    = dn.*D;
    Dn{n} = dn.*J;
end


% initialise (solve matrix polynomial for Q = QJ)
% Q = sum( ((D.^n).*J)*(Q^n)/factorial(n) ): n = 0,1,...
%==========================================================================
for i = 1:N
    
    QJn   = 1;
    Q     = J;
    for n = 1:N
        
        % n-th order Taylor term
        %------------------------------------------------------------------
        QJn = QJn*QJ/n;
        dQ  = Dn{n}*QJn;
        Q   = Q + dQ;
        
        % break if convergence
        %------------------------------------------------------------------
        if norm(dQ,'inf') < TOL; break, end
        
    end
    
    % break if unstable
    %----------------------------------------------------------------------
    if any(any(isnan(Q))), Q = QJ; break, end
    
    % Robbins–Monro update and break if convergence
    %----------------------------------------------------------------------
    QJ = QJ*(1 - a) + Q*a;
    
    if norm((QJ - Q),'inf') < TOL; break, end
    
end
Q      = Q*spm_inv(J);
