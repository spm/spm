function X = spm_inv(A,TOL)
% Inverse for ill-conditioned matrices
% FORMAT X = spm_inv(A,TOL)
%
% A   - matrix
% X   - inverse
%
% TOL - tolerance: default = max(eps(norm(A,'inf'))*max(m,n),exp(-32))
%
% This routine simply adds a small diagonal matrix to A and calls inv.m
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end
 
% tolerance
%--------------------------------------------------------------------------
if nargin == 1
    TOL  = max(eps(norm(A,'inf'))*max(m,n),exp(-32));
end

% inverse
%--------------------------------------------------------------------------
X     = inv(A + speye(m,n)*TOL);
