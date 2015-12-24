function [Y] = spm_dot(X,x,DIM)
% Multidimensional dot (inner) preoduct
% FORMAT [Y] = spm_dot(X,x,DIM)
%
% X  - numeric array
% x  - vector or cell array of numeric vectors
%
% Y  - inner product obtained by summing the products of X and x along DIM
%
% If DIM is not specified the last dimension of X is used.  If x is a cell
% array recursive dot products are computed (starting with the last entry
% if (the vector) DIM is not specified).
%
% See also: spm_cross
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dot.m 6656 2015-12-24 16:49:52Z guillaume $


% initialise X and vX
%--------------------------------------------------------------------------
if iscell(x)
    if nargin < 3
        DIM = (1:numel(x)) + numel(size(X)) - numel(x);
    end
    for i = 1:numel(x)
        X   = spm_dot(X,x{i},DIM(i));
        DIM = DIM - 1;
    end
    Y     = X;
    return
end

% inner product
%==========================================================================
if nargin < 3, DIM = numel(size(X)); end

% deal with simple cases
%--------------------------------------------------------------------------
if isvector(X)
    Y = X*x(:);
    return
elseif ismatrix(X)
    if DIM == 1
        Y = x'*X;
    else
        Y = X*x;
    end
    return
end

d      = size(X);
ind    = cell(size(d));
ind(:) = {':'};
Y      = X;

for i = 1:numel(x)
    sub       = ind;
    sub{DIM}  = i;
    Y(sub{:}) = X(sub{:}).*x(i);
end
Y = sum(Y,DIM);
if ~ismatrix(Y)
    d(DIM) = [];
    Y      = reshape(Y,d);
end
