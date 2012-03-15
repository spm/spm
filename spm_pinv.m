function X = spm_pinv(A)
% pseudoinverse for sparse matrices
% FORMAT X = spm_pinv(A)
%
% X   - matrix
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_pinv.m 4690 2012-03-15 20:10:40Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end


% try generalised inverse
%--------------------------------------------------------------------------
sw = warning('off','MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
X     = spm_inv(A'*A);
warning(sw);

if all(isfinite(X(:)))
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

