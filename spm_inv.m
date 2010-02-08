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
% $Id: spm_inv.m 3716 2010-02-08 13:58:09Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end
 
% inverse
%--------------------------------------------------------------------------
TOL = eps(norm(A,1))*max(m,n);
X   = inv(A + speye(m,n)*TOL);
