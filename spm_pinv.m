function X = spm_pinv(A,TOL)
% pseudo-inverse for sparse matrices
% FORMAT X = spm_pinv(A,TOL)
%
% A   - matrix
% TOL - Tolerance to force singular value decomposition
% X   - generalised inverse
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_pinv.m 5877 2014-02-11 20:03:34Z karl $
 
% check A
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end
 
 
% try generalised inverse
%--------------------------------------------------------------------------
if nargin < 2
    sw = warning('off','MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    X  = spm_inv(A'*A);
    warning(sw);
    
    
    % check everything is finite
    %----------------------------------------------------------------------
    if all(isfinite(X(:)))
        X = X*A';
        return
    end
end
 
% pseudo-inverse
%--------------------------------------------------------------------------
[U,S,V] = spm_svd(A,0);
S       = full(diag(S));

if nargin < 2, TOL = max(m,n)*eps(max(S)); end

% tolerance
%------------------------------------------------------_-------------------
r   = sum(abs(S) > TOL);
if ~r
    X = sparse(n,m);
else
    i = 1:r;
    S = sparse(i,i,1./S(i),r,r);
    X = V(:,i)*S*U(:,i)';
end
