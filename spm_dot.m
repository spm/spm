function [X] = spm_dot(X,x,i)
% Multidimensional dot (inner) preoduct
% FORMAT [Y] = spm_dot(X,x,[DIM])
%
% X   - numeric array
% x   - cell array of numeric vectors
% DIM - dimensions to omit (asumes ndims(X) = numel(x))
%
% Y  - inner product obtained by summing the products of X and x along DIM
%
% If DIM is not specified the leading dimensions of X are omitted.
%
% See also: spm_cross
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dot.m 6856 2016-08-10 17:55:05Z karl $

% initialise X and vXthere
%--------------------------------------------------------------------------
if nargin < 3
    DIM    = (1:numel(x)) + ndims(X) - numel(x);
else
    DIM    = (1:numel(x)) + ndims(X) - numel(x);
    DIM(i) = [];
    x(i)   = [];
end

% inner product using bsxfun
%----------------------------------------------------------------------
for d = 1:numel(x)
    s         = ones(1,ndims(X));
    s(DIM(d)) = numel(x{d});
    X         = bsxfun(@times,X,reshape(full(x{d}),s));
    X         = sum(X,DIM(d));
end

% eliminate Singleton dimensions
%----------------------------------------------------------------------
X = squeeze(X);



