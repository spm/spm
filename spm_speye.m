function [D] = spm_speye(m,n,k,c)
% Sparse leading diagonal matrix
% FORMAT [D] = spm_speye(m,n,k,c)
%
% returns an m x n matrix with ones along the k-th leading diagonal. If
% called with an optional fourth argument c = 1, a wraparound sparse matrix
% is returned. If c = 2, then empty rows or columns are filled in on the
% leading diagonal.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging


% default k = 0
%--------------------------------------------------------------------------
if nargin < 4, c = 0; end
if nargin < 3, k = 0; end
if nargin < 2, n = m; end
 
% leading diagonal matrix
%--------------------------------------------------------------------------
D = spdiags(ones(m,1),k,m,n);

% add wraparound if necessary
%--------------------------------------------------------------------------
if c == 1
    if k < 0
        D = D + spm_speye(m,n,min(n,m) + k);
    elseif k > 0
        D = D + spm_speye(m,n,k - min(n,m));
    end
elseif c == 2
    i = find(~any(D));
    D = D + sparse(i,i,1,n,m);
    
end
