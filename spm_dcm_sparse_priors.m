function [A,K,k] = spm_dcm_sparse_priors(n)
% Return Adjacency matrices for bidirectional coupling
% FORMAT  [A,K,k] = spm_dcm_sparse_priors(n)
%
% INPUT:
%    n         - number of nodes
%
% OUTPUT:
%    A{:}      - adjacency matrices
%    K{1:K}{:} - adjacency matrices (for k - 1 edges)
%    k         - row vector of edge numbers (size)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
if nargin < 2; k = n - 1; end

% get permutations of connections
%==========================================================================
N     = (n*n - n)/2;                % number of (bidirectional) connections
np    = [1; 0];
for i = 1:N
    p     = np;
    np    = [p; p];
    for j = 1:i
       np(:,j)    = kron(p(:,j),[1; 1]);
    end
    np(:,(i + 1)) = kron(ones(2^i,1),[1; 0]);
end


% get permutations of K connections (in A{:})
%==========================================================================
[I,J] = find(triu(ones(n,n),1));
K     = cell(N + 1,1);
A     = {};
for i = 1:2^N
    
    % adjacency matrix
    %----------------------------------------------------------------------
    a = full(sparse(I,J,p(i,:),n,n));
    a = a + a' + eye(n,n);
    
    % save in arrays
    %----------------------------------------------------------------------
    k(i,1)               = sum(p(i,:));
    K{1 + k(i)}{end + 1} = a;
    A{i}                 = a;
    
end
