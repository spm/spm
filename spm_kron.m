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
%   Class support for inputs X,Y:
%      float: double, single
% 
% When called with a single cell array input, the tensor product
% is formed recursively 

%   Previous versions by Paul Fackler, North Carolina State,
%   and Jordan Rosenthal, Georgia Tech.
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 7811 $ $Date: 2004/06/25 18:52:18 $
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_kron.m 7811 2020-04-05 12:00:43Z karl $


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
K = matlab.internal.sparse.kronSparse(sparse(A), sparse(B));


