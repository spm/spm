function [M0,q0,X,x,f,M1,L] = spm_fp(M,x,u)
% Fokker-Planck operators and equilibrium density for dynamic systems
% FORMAT [M0,q0,X,x,f,M1,L] = spm_fp(M,x,u)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt    = f(x,u,P) or f(x,u,a,P)  {function string or m-file}
%    M.g   - y(t)     = g(x,u,P)                {function string or m-file}
%    M.m   - m inputs
%    M.n   - n states
%    M.l   - l outputs
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - precision matrix of state noise
% x    - cell array of vectors specifying evaluation grid
% u    - expansion point for inputs or causes;
%
% M0   - 1st order FP operator dq/dt = M0*q + u*M1*q,  q = p(X);
% q0   - stable or equilibrium mode: M0*q0 = 0
% X    - evaluation points of state space
% x    - cell array of vectors specifying evaluation grid
% f    - flow
% M1   - 2nd order FP operator
% L    - output matrix                 <y> = L*q;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp.m 5219 2013-01-29 17:07:07Z spm $
 
% default: first level of hierarchical model
%--------------------------------------------------------------------------
M   = M(1);
 
% default expansion point for inputs
%--------------------------------------------------------------------------
if nargin < 3
    u0 = zeros(M.m,1);
else
    u0 = u;
end

% event space: get or create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
n  = length(M.x);
if nargin < 2

    % use M.X
    %----------------------------------------------------------------------
    X     = M.X;
    for i = 1:size(X,2)
        x{i} = unique(X(:,i));
    end
else
    
    % use x
    %----------------------------------------------------------------------
    [X,x] = spm_ndgrid(x);
end

% f(x,0)
%--------------------------------------------------------------------------
N     = length(X);
f     = sparse(n,N);
try
    for i = 1:N
        f(:,i) = feval(M.f,X(i,:)',u0,M.pE);
    end
catch
    for i = 1:N
        f(:,i) = feval(M.f,X(i,:)',u0,0,M.pE);
    end
end


% Fokker-Planck operator; i.e., Jacobian J - M0
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
        Dp = spm_kron(I{j},Dp);
        Dq = spm_kron(I{j},Dq);
    end
    Dp    = spm_kron(dp,Dp);
    Dq    = spm_kron(dq,Dq);
    for j = (i + 1):n
        Dp = spm_kron(I{j},Dp);
        Dq = spm_kron(I{j},Dq);
    end
    DP{i} = Dp;
 
    % augment Jacobian
    %----------------------------------------------------------------------
    J     = J + Dp*spdiags(u,0,N,N) + Dq*spdiags(v,0,N,N);
 
end
 
% dispersion
%--------------------------------------------------------------------------
C     = inv(M.W);
if length(C) ~= n
    C = C(1)*speye(n,n);
end
for i = 1:n
    for j = 1:n
        if C(i,j)
            J  = J - C(i,j)*DP{i}*DP{j}'/2;
        end
    end
end
 
% return if only M0 is required
%--------------------------------------------------------------------------
M0    = J;
if nargout == 1, return, end

 
% stable mode - stochastic iteration
%==========================================================================
q     = sparse(N,1) + 1/N;
for i = 1:1024
    dq = M0*q;
    if norm(dq,1) < exp(-4), break, end
    while min(q + dq) < 0
        dq = dq - dq/2;
    end
    q  = q + dq;
end
for i = 1:n
    nx{i} = length(x{i});
end
q0    = reshape(q,nx{:});

if nargout < 6,  return,  end

% input induced changes to Jacobian dJ/du - M1
%==========================================================================
M1    = spm_diff('spm_fp',M,x,u0,3);
if nargout < 6,  return,  end
 
 
% Output matrix L - M1
%==========================================================================
 
% l(x,0)
%--------------------------------------------------------------------------
L     = sparse(M.l,N);
for i = 1:N
    try
        L(:,i) = feval(M.g,X(i,:)',u0,M.pE);
    catch
        L(:,i) = feval(M.g,X(i,:)',u0,[],M.pE);
    end
end

return

% NOTES: analytic solution for q0
%--------------------------------------------------------------------------
q   = length(M0);
J   = [M0; ones(1,q)];
q0  = sum(inv(J'*J),2);
