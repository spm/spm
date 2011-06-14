function X = spm_inv_spd(A, TOL)
% inverse for symmetric positive (semi)definite matrices
% FORMAT X = spm_inv_spd(A,TOL)
%
% A   - symmetric positive definite matrix (e.g. covariance or precision)
% X   - inverse (should remain symmetric positive definite)
%
% TOL - tolerance: default = exp(-32)
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Ged Ridgway
% $Id: spm_inv_spd.m 4360 2011-06-14 16:46:37Z ged $

% if ~all(isfinite(A(:))), error('Matrix has non-finite elements!'); end

if nargin < 2
    TOL  = exp(-32);
end

[i j] = find(A);
if isempty(i)
    % Special cases:  empty or all-zero matrix, return identity/TOL
    %----------------------------------------------------------------------
    X = eye(length(A)) / TOL;
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
elseif norm(A - A', 1) < TOL
    % symmetric, try LDL factorisation (but with L->X to save memory)
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
    else
        error('Matrix is not positive semi-definite according to ldl')
    end
else
    error('Matrix is not symmetric to given tolerance');
end

% if ~all(isfinite(X(:))), error('Inverse has non-finite elements!'); end

function d = invtol(d, TOL)
% compute reciprocal of values, clamped to lie between TOL and 1/TOL
if any(d < -TOL)
    error('Matrix is not positive semi-definite at given tolerance')
end
d = max(d, TOL);
d = 1./d;
d = max(d, TOL);
