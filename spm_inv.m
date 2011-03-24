function X = spm_inv(A)
% inverse for ill-conditioned matrices
% FORMAT X = spm_inv(A)
%
% A  - matrix
% X  - inverse
%
% This routine simply adds a small identity matrix to A and calls inv.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_inv.m 4261 2011-03-24 16:39:42Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end
 
% tolerance
%--------------------------------------------------------------------------
TOL = max(eps(norm(A,'inf'))*max(m,n),exp(-32));

% inverse
%--------------------------------------------------------------------------
X   = inv(sparse(A) + speye(m,n)*TOL);
