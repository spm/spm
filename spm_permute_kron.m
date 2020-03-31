function A = spm_permute_kron(A,dim,order)
% permutation of a Kronecker tensor product 
% FORMAT A = spm_kron(A,DIM,ORDER)
% A     - n-dimensional array (A1 x A2 x ...
% dim   - dimensions [length(A1), length(A2), ...
% order - re-ordering; e.g., [2,1, ...
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_permute_kron.m 7809 2020-03-31 11:55:09Z karl $


% creating indices
%--------------------------------------------------------------------------
i     = reshape(1:prod(dim),dim);
i     = permute(i,order);
A     = A(i,i);

