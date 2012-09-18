function [M0,M1,L,X,q0] = spm_mfa(M,x,u)
% Jacobian for mean field approximations
% FORMAT [M0,M1,L,X,q0] = spm_mfa(M,x,u)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt    = f(x,u,P)                     {function string or m-file}
%    M.g   - y(t)     = g(x,u,P)                     {function string or m-file}
%    M.m   - m inputs
%    M.n   - n states
%    M.l   - l ouputs
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - covariance matrix of deterministic noise
% x    - cell array of vectors specifying evaluation grid
% u    - expansion point for inputs (c.f. background activity);
%
% M0   - 1st order Bilinear matrix dq/dt = M0*q + u*M1*q,  q = p(X);
% M1   - 2nd order Bilinear matrix
% L    - output matrix                 <y> = L*q;
% X    - evaluation points of state space
% q0   - stable mode M0*q0 = 0
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfa.m 4936 2012-09-18 19:47:55Z karl $
 
 
 
% default expansion point for inputs
%--------------------------------------------------------------------------
if nargin < 3
    u0 = zeros(M.m,1);
else
    u0 = u;
end
 
% create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
n     = length(x);
for i = 1:n
    q     = 1;
    for j = 1:n
        q = kron(x{j}(:).^(i == j),q);
    end
    X(:,i) = q;
end
 
% f(x,0)
%--------------------------------------------------------------------------
N     = length(X);
f     = sparse(n,N);
for i = 1:N
    f(:,i) = feval(M.f,X(i,:)',u0,M.pE);
end
 
% Jacobian J - M0
%==========================================================================
for i = 1:n
    d(i)  = length(x{i});
    dx(i) = x{i}(2) - x{i}(1);
    I{i}  = speye(d(i),d(i));
end
J     = sparse(N,N);
 
 
% cycle over dimensions of state space
%--------------------------------------------------------------------------
for i = 1:length(x)
 
    % differential operators (for positive and negative flows)
    %----------------------------------------------------------------------
    j          = find(f(i,:) < 0);
    u          = sparse(j,1,f(i,j),N,1);
    j          = find(f(i,:) > 0);
    v          = sparse(j,1,f(i,j),N,1);
 
    dp         = spdiags(ones(d(i),1)*[1 -1],[ 0;1],d(i),d(i))/dx(i);
    dp(:,1)    = 0;
    dq         = spdiags(ones(d(i),1)*[1 -1],[-1;0],d(i),d(i))/dx(i);
    dq(:,end)  = 0;
 
    % Kronecker tensor products
    %----------------------------------------------------------------------
    Dp         = 1;
    Dq         = 1;
    for j = 1:(i - 1)
        Dp = kron(I{j},Dp);
        Dq = kron(I{j},Dq);
    end
    Dp    = kron(dp,Dp);
    Dq    = kron(dq,Dq);
    for j = (i + 1):n
        Dp = kron(I{j},Dp);
        Dq = kron(I{j},Dq);
    end
    DP{i} = Dp;
 
    % augment Jacobian
    %----------------------------------------------------------------------
    J     = J + Dp*spdiags(u,0,N,N) + Dq*spdiags(v,0,N,N);
 
end
 
% dispersion
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if M.W(i,j)
            J  = J - M.W(i,j)*DP{i}*DP{j}'/2;
        end
 
    end
end
 
% return if only M0 is required
%--------------------------------------------------------------------------
M0  = J;
if nargout == 1,  return,  end
 
 
% input induced changes to Jacobian dJ/du - M1
%==========================================================================
du    = 1e-6;
M1    = {};
for i = 1:M.m
    u     = u0;
    u(i)  = du;
    Mu    = spm_mfa(M,x,u);
    M1{i} = (Mu - M0)/du;
end
 
 
% return if only M0, M1 are required
%--------------------------------------------------------------------------
if nargout == 2,  return,  end
 
 
% Output matrix L - M1
%==========================================================================
 
% l(x,0)
%--------------------------------------------------------------------------
L     = sparse(M.l,N);
for i = 1:N
    L(:,i) = feval(M.g,X(i,:)',u,M.pE);
end
 
% equilibrium mode
%--------------------------------------------------------------------------
if nargout == 5
 
    % stable mode - stochastic iteration
    %----------------------------------------------------------------------
    q   = ones(N,1)/N;
    for i = 1:1024
        dq = M0*q;
        q  = q + dq/N;
        q  = q/sum(q);
    end
    q0  = q;
 
    % analytic solution (not implemented)
    %----------------------------------------------------------------------
    % q   = length(M0);
    % J   = [M0; ones(1,q)];
    % q0  = sum(inv(J'*J),2);
 
end
