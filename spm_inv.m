function X = spm_inv(A,TOL)
% inverse for ill-conditioned matrices
% FORMAT X = spm_inv(A,TOL)
%
% A   - matrix
% X   - inverse
%
% TOL - tolerance: default = max(eps(norm(A,'inf'))*max(m,n),exp(-32))
%
% This routine simply adds a small diagonal matrix to A and calls inv.m
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston and Ged Ridgway
% $Id: spm_inv.m 4317 2011-04-27 20:39:03Z ged $
 
if nargin < 2
    TOL  = max(eps(norm(A,'inf'))*length(A), exp(-32));
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
    X = diag(d);
elseif norm(A - A', 1) < TOL
    % symmetric, use LDL factorisation (but with L->X to save memory)
    %----------------------------------------------------------------------
    [X D P] = ldl(full(A)); % P'*A*P = L*D*L', A = P*L*D*L'*P'
    d = diag(D);
    d = invtol(d, TOL);
    % inv(A) = P*inv(L')*inv(D)*inv(L)*P' = (L\P')'*inv(D)*(L\P')
    % triangular system should be quick to solve and stay approx triangular
    X = X\P';
    X = X'*diag(d)*X;
else
    % asymmetric matrix -- use old approach (SVD seems not to work well)
    %----------------------------------------------------------------------
    X = inv(A + TOL*eye(length(A)));
    % [U S X] = svd(full(A));
    % d = diag(S);
    % d = invtol(d, TOL);
    % X = X*diag(d)*U';
end
if issparse(A), X = sparse(X); end
    
function d = invtol(d, TOL)
% compute reciprocal of values, clamped to lie between TOL and 1/TOL
d = max(d, TOL);
d = 1./d;
d = max(d, TOL);
