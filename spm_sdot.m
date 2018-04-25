function [Y] = spm_sdot(X,x,k)
% Sparse multidimensional dot (inner) product
% FORMAT [Y] = spm_sdot(X,x,[DIM])
%
% X   - numeric array
% x   - cell array of numeric vectors
% DIM - dimension to omit (asumes ndims(X) = numel(x))
%
% Y  - inner product obtained by summing the products of X and x along DIM
%
% If DIM is not specified the leading dimensions of X are omitted. This
% routine assumes X is sparse
%
% See also: spm_dot, spm_cross
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sdot.m 7300 2018-04-25 21:14:07Z karl $

% initialise dimensions
%--------------------------------------------------------------------------
DIMS   = size(X);
NDIM   = numel(DIMS);
J      = 1:NDIM;
J(k)   = [];
I      = cell(1,NDIM);
j      = find(X > exp(-16));
[I{:}] = ind2sub(DIMS,j);

% sum of products
%--------------------------------------------------------------------------
Y      = zeros(numel(x{k}),1);
for i  = 1:numel(j)
    p  = X(j(i));
    for d = J
        p = p*x{d}(I{d}(i));
    end
    Y(I{k}(i)) = Y(I{k}(i)) + p;
end



