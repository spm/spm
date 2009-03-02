function X = spm_pinv(A)
% psudoinverse for sparse matrices
% FORMAT X = spm_pinv(A)
%
% X   - matrix
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_pinv.m 2804 2009-03-02 12:03:00Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end


% try generalised inverse
%--------------------------------------------------------------------------
warning off
X     = inv(A'*A);
warning on

if ~any(isnan(diag(X)))
    X = X*A';
    return
end

% pseudoinverse
%--------------------------------------------------------------------------
[U,S,V] = spm_svd(A,0);

S   = full(diag(S));
TOL = max(m,n)*eps(max(S));
r   = sum(abs(S) > TOL);
if ~r
    X = sparse(n,m);
else
    i = 1:r;
    S = sparse(i,i,1./S(i),r,r);
    X = V(:,i)*S*U(:,i)';
end

