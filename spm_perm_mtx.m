function [K] = spm_perm_mtx(n)
% Return a matrix of indices permuted over n
% FORMAT [K] = spm_perm_mtx(n)
%    n   - (scalar) number of indices
%    K   - (2^n x n) permutation matrix
% or
%    n   - (vector) indices
%    K   - (length(n)! x n) permutation matrix
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

 
% get permutations
%==========================================================================

% permute zeros and ones
%--------------------------------------------------------------------------
if isscalar(n)
    
    N  = 2^n;
    K  = sparse(N,n);
    x  = sparse(1,1,1,2,1);
    for i = 1:n
        y      = ones(N/length(x),1);
        K(:,i) = kron(x,y);
        x      = [x;x];
    end
    
% permute indices
%--------------------------------------------------------------------------
elseif isvector(n)
    
    n  = n(:);
    K  = n;
    while size(K,2) < length(n)
        x     = [];
        for i = 1:size(K,1)
            d = K(i,:);
            r = n;
            for j = 1:length(d)
                r(r == d(j)) = [];
            end
            x = [x; [kron(ones(length(r),1),d) r]];
        end
        K = x;
    end
end

% make logical
%--------------------------------------------------------------------------
K  = logical(K);
