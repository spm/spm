function [K] = spm_perm_mtx(n)
% Returns a matrix of indices permuted over n
% FORMAT [K] = spm_perm_mtx(n)
%    n   - (scalar) number of indices
%    K   - (n^2 x n) permutation matrix
%    n   - (vector) indices
%    K   - (length(n)! x n) permutation matrix
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_perm_mtx.m 4100 2010-10-22 19:49:17Z karl $
 
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