function A = spm_permute_kron(A,dim,order)
% Permutation of a Kronecker tensor product 
% FORMAT A = spm_permute_kron(A,dim,order)
% A     - 2-dimensional array (A1 x A2 x ...
% dim   - dimensions [length(A1), length(A2), ...
% order - re-ordering; e.g., [2,1, ...
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


% create indices
%--------------------------------------------------------------------------
i     = reshape(1:prod(dim),dim);
i     = permute(i,order);
A     = A(i,i);
