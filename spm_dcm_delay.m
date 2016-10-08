function [Q,J] = spm_dcm_delay(P,M,J,N)
% returns the delay operator for flow and Jacobians of dynamical systems
% FORMAT [Q,J] = spm_dcm_delay(P,M)
% P   - model parameters
% M   - model specification structure
% J   - optional: system Jacobian
% N   - optional: auto Taylor expansion [default: 2^6]
%
% Required fields:
%   M.f - dx/dt    = f(x,u,P,M)            {function string or m-file}
%   M.m - m inputs
%   M.n - n states
%   M.x - (n x 1) = x(0) = expansion point: defaults to x = 0;
%   M.u - (m x 1) = u    = expansion point: defaults to u = 0;
%
%
% Returns the delay operator for Jacobians of dynamical systems where the
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
% $Id: spm_dcm_delay.m 6900 2016-10-08 13:16:46Z karl $

% order of Taylor approximation
%--------------------------------------------------------------------------
if nargin < 4, N = 2^8; end

% Jacobian
%==========================================================================
if nargin < 3
    
    if isfield(M,'u'), u = spm_vec(M.u); else, u = sparse(M.m,1); end
    
    J = full(spm_diff(M.f,M.x,u,P,M,1));
    
else
    J = full(J);
end

% evaluate delay matrix D from parameters
%==========================================================================

% paramterised delays
%--------------------------------------------------------------------------
if isfield(P,'D')
    
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
    
    if isnumeric(M.x)
        
        % number of states per sources
        %------------------------------------------------------------------
        nx  = size(M.x,2);
        De  = kron(ones(nx,nx),De);
        Di  = kron(ones(nx,nx) - speye(nx,nx),Di);
        D   = Di + De;
        
    else
        
        % number of states per sources
        %------------------------------------------------------------------
        for i = 1:numel(M.x)
            for j = 1:numel(M.x)
                m = numel(M.x{i});
                n = numel(M.x{j});
                if i == j
                    D{i,j} = Di(i,j)*(1 - speye(m,n));
                else
                    D{i,j} = De(i,j)*(1 - zeros(m,n));
                end
            end
        end
        D   = spm_cat(D);
        
    end
   
else
    
    % no delays
    %----------------------------------------------------------------------
    D = sparse(0);
    
end

% suppress delays between voltage and current
%--------------------------------------------------------------------------
if isfield(M,'nodelay'), D(J == 1) = 0; end


% Jacobian and delay operator
%==========================================================================

% delay operator: first-order approximation if N = 0
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
if ~N
    Q  = inv(speye(length(J)) + D.*J);
    return
end

% delay operator: estimated using a Robbins–Monro algorithm
%--------------------------------------------------------------------------
D     = -D;
QJ    = (eye(length(J)) - D.*J)\J;
a     = 1/4;
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
        QJn           = QJn*QJ/n;
        dQ            = Dn{n}*QJn;
        dQ(isnan(dQ)) = 0;
        Q             = Q + dQ;
        
        % break if convergence
        %------------------------------------------------------------------
        if norm(dQ,'inf') < TOL; break, end
        
    end

    % Robbins–Monro update and break if convergence
    %----------------------------------------------------------------------
    QJ = QJ*(1 - a) + Q*a;
    
    if norm((QJ - Q),'inf') < TOL; break, end
    
end
Q      = Q*spm_inv(J);
