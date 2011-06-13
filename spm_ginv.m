function X = spm_ginv(A,TOL)
% inverse for ill-conditioned matrices
% FORMAT X = spm_ginv(A,TOL)
%
% A   - matrix
% X   - inverse
%
% TOL - tolerance: default = exp(-32)
%
% This routine simply adds a small diagonal matrix to A and calls inv.m
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston and Ged Ridgway
% $Id: spm_ginv.m 4354 2011-06-13 18:53:11Z karl $

% if ~all(isfinite(A(:))), error('Matrix has non-finite elements!'); end

if nargin < 2
    TOL  = exp(-32);
end

[i j] = find(A);
if isempty(i)
    % Special cases:  empty or all-zero matrix, return identity/TOL
    %----------------------------------------------------------------------
    X = eye(length(A)) / TOL;
    return
elseif all(i == j)
    % diagonal matrix
    %----------------------------------------------------------------------
    d = diag(A);
    d = invtol(d, TOL);
    if issparse(A)
        n = length(A);
        X = sparse(1:n, 1:n, d);
    else
        X = diag(d);
    end
    return
elseif norm(A - A', 1) < TOL
    % symmetric, use LDL factorisation (but with L->X to save memory)
    %----------------------------------------------------------------------
    [X D P] = ldl(full(A)); % P'*A*P = L*D*L', A = P*L*D*L'*P'
    [i j d] = find(D);
    % non-diagonal values indicate not positive semi-definite
    if all(i == j)
        d = invtol(d, TOL);
        % inv(A) = P*inv(L')*inv(D)*inv(L)*P' = (L\P')'*inv(D)*(L\P')
        % triangular system should be quick to solve and stay approx tri.
        X = X\P';
        X = X'*diag(d)*X;
        if issparse(A), X = sparse(X); end
        return
    end
end

% If still here, either asymmetric or non-diagonal and not pos. semi-def.
%--------------------------------------------------------------------------
% Approach from original spm_inv (not ideal for negative eigenvalues...)
X = inv(A + TOL*speye(length(A))); % returns sparse if A was, else full
% SVD-based approach (much slower, and seemingly not much better...)
% [U S X] = svd(full(A));
% d = diag(S);
% d = invtol(d, TOL);
% X = X*diag(d)*U';
% if issparse(A), X = sparse(X); end

% if ~all(isfinite(X(:))), error('Inverse has non-finite elements!'); end

function d = invtol(d, TOL)
% compute reciprocal of values, clamped to lie between TOL and 1/TOL,
% although if any values below -TOL, keep their negative signs.
s = d < -TOL;
d = max(abs(d), TOL);
d = 1./d;
d = max(d, TOL);
d(s) = -d(s);