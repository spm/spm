function [U,S,V] = spm_svd(X,U,T)
% computationally efficient SVD (that can handle sparse arguments)
% FORMAT [U,S,V] = spm_svd(X,u,t);
% X    - {m x n} matrix
% u    - threshold for normalized eigenvalues (default = 1e-6)
% t    - threshold for raw eigenvalues        (default = 0)
%
% U    - {m x p} singular vectors
% V    - {m x p} singular variates
% S    - {p x p} singular values
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_svd.m 4124 2010-11-18 16:56:53Z karl $



% default thresholds
%---------------------------------------------------------------------------
if nargin < 2, U = 1e-6; end
if nargin < 3, T = 0;    end

% deal with sparse matrices
%---------------------------------------------------------------------------
[M N] = size(X);
p     = find(any(X,2));
q     = find(any(X,1));
X     = X(p,q);

% SVD
%---------------------------------------------------------------------------
[i j s] = find(X);
[m n]   = size(X);
if any(i - j)

    % off-leading diagonal elements - full SVD
    %-------------------------------------------------------------------
    X     = full(X);
    if m > n

        [v S v] = svd(X'*X,0);
        S       = sparse(S);
        s       = diag(S);
        j       = find(s*length(s)/sum(s) >= U & s >= T);
        v       = v(:,j);
        u       = spm_en(X*v);
        S       = sqrt(S(j,j));

    elseif m < n

        [u S u] = svd(X*X',0);
        S       = sparse(S);
        s       = diag(S);
        j       = find(s*length(s)/sum(s) >= U & s >= T);
        u       = u(:,j);
        v       = spm_en(X'*u);
        S       = sqrt(S(j,j));

    else

        [u S v] = svd(X,0);
        S       = sparse(S);
        s       = diag(S).^2;
        j       = find(s*length(s)/sum(s) >= U & s >= T);
        v       = v(:,j);
        u       = u(:,j);
        S       = S(j,j);
    end

else
    S     = sparse(1:n,1:n,s,m,n);
    u     = speye(m,n);
    v     = speye(m,n);
    [i j] = sort(-s);
    S     = S(j,j);
    v     = v(:,j);
    u     = u(:,j);
    s     = diag(S).^2;
    j     = find(s*length(s)/sum(s) >= U & s >= T);
    v     = v(:,j);
    u     = u(:,j);
    S     = S(j,j);

end

% replace in full matrices
%---------------------------------------------------------------------------
j      = length(j);
U      = sparse(M,j);
V      = sparse(N,j);
if j
    U(p,:) = u;
    V(q,:) = v;
end
