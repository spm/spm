function [X] = spm_dot(X,x,i)
% Multidimensional dot (inner) product
% FORMAT [Y] = spm_dot(X,x,[DIM])
%
% X   - numeric array
% x   - cell array of numeric vectors
% DIM - dimensions to omit (assumes ndims(X) = numel(x))
%
% Y  - inner product obtained by summing the products of X and x along DIM
%
% If DIM is not specified the leading dimensions of X are omitted.
% If x is a vector the inner product is over the first matching dimension of X
%
% See also: spm_cross
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging

% scalar product
%--------------------------------------------------------------------------
if numel(spm_vec(x)) == 1
    X = X*spm_vec(x);
    return
end

% initialise dimensions
%--------------------------------------------------------------------------
if iscell(x)

    % omit leading dimensions
    %----------------------------------------------------------------------
    DIM = (1:numel(x)) + max(ndims(X),numel(x)) - numel(x);
    
else

    %  find first matching dimension
    %----------------------------------------------------------------------
    DIM = find(size(X) == numel(x),1);
    x   = {x};

end

% omit dimensions specified
%--------------------------------------------------------------------------
if nargin > 2
    DIM(i) = [];
    x(i)   = [];
end

% inner product using tensorprod
%--------------------------------------------------------------------------
for d = 1:numel(x)
    X    = tensorprod(X,full(x{d}),DIM(d),1);
    DIM  = DIM - 1;
end

return

