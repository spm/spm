function K = spm_kron(A,B)
% Kronecker tensor product with sparse outputs 
% FORMAT K = spm_kron(A,B)
%        K = spm_kron(A)
%
%   KRON(X,Y) is the Kronecker tensor product of X and Y.
%   The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y.   For
%   example, if X is 2 by 3, then KRON(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
% 
% When called with a single cell array input, the tensor product
% is formed recursively 
%__________________________________________________________________________
% Copyright (C) 2008-2020 Wellcome Centre for Human Neuroimaging
 
% Karl Friston
% $Id: spm_kron.m 7814 2020-04-06 11:47:47Z guillaume $


% Deal with cell arrays
%--------------------------------------------------------------------------
if iscell(A)
    K = 1;
    for i = 1:numel(A)
        K = kron(A{i},K);
    end
    return
end

% Kronecker tensor product
%--------------------------------------------------------------------------
K = kron(sparse(A), sparse(B));
