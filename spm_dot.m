function [X] = spm_dot(X,x,DIM)
% Multidimensional dot (inner) preoduct
% FORMAT [Y] = spm_dot(X,x,DIM)
%
% X  - numeric array
% x  - cell array of numeric vectors
%
% Y  - inner product obtained by summing the products of X and x along DIM
%
% If DIM is not specified the trailing dimensions of X are used.
%
% See also: spm_cross
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dot.m 6812 2016-06-18 11:16:21Z karl $

% initialise X and vXthere
%--------------------------------------------------------------------------
if nargin < 3
    DIM = (1:numel(x)) + numel(size(X)) - numel(x);
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



