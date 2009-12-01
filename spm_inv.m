function X = spm_inv(A)
% inverse for ill-conditioned matrices
% FORMAT X = spm_inv(A)
%
% X   - matrix
%
% This routine simply adds a small identity matrix to A and calls inv.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_inv.m 3605 2009-12-01 13:29:43Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end
 
% inverse
%--------------------------------------------------------------------------
TOL = max(m,n)*eps(norm(A,1));
X   = inv(A + speye(m,n)*TOL);
