function A = spm_permute_kron(A,dim,order)
% permutation of a Kronecker tensor product 
% FORMAT A = spm_permute_kron(A,dim,order)
% A     - 2-dimensional array (A1 x A2 x ...
% dim   - dimensions [length(A1), length(A2), ...
% order - re-ordering; e.g., [2,1, ...
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_permute_kron.m 8171 2021-10-25 10:14:50Z karl $


% create indices
%--------------------------------------------------------------------------
i     = reshape(1:prod(dim),dim);
i     = permute(i,order);
A     = A(i,i);

