function [X] = spm_dot(X,x,i)
% Multidimensional dot (inner) product
% FORMAT [Y] = spm_dot(X,x,[DIM])
%
% X   - numeric array
% x   - numeric vector or cell array of vectors
% DIM - dimensions to skip [assumes ndims(X) = numel(x)]
%
% Y   - inner product obtained by summing the products of X and x
%
% If DIM is not specified the leading dimensions of X are skipped. If x is
% a vector the inner product is over the first (matching) dimension of X.
% This means that if called with a vector valued x, the dot product will be
% over the first (matching) dimension. Conversely, if called with {x} the
% dot product will be over the last (matching) dimension of X.
%
% This version calls tensorprod.m
%
% See also: spm_cross
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging


% initialise dimensions
%--------------------------------------------------------------------------
if iscell(x)

    % scalar product
    %--------------------------------------------------------------------------
    if numel(x) == 1
        if numel(x{1}) == 1
            X = X*x{1};
            return
        end
    end

    % omit leading dimensions
    %----------------------------------------------------------------------
    DIM = (1:numel(x)) + max(ndims(X),numel(x)) - numel(x);
    
else

    % scalar product
    %--------------------------------------------------------------------------
    if numel(x) == 1
        X = X*x;
        return
    end

    % find first matching dimension
    %----------------------------------------------------------------------
    DIM = find(size(X) == numel(x),1);
    x   = {x};

end

% omit specified dimensions 
%--------------------------------------------------------------------------
if nargin > 2
    DIM(i) = [];
    x(i)   = [];
end

% inner product using tensorprod
%--------------------------------------------------------------------------
if islogical(X),X = single(X); end
if issparse(X), X = full(X);   end
for d = 1:numel(x)
    
    if islogical(x{d}),x{d} = single(x{d}); end
    if issparse(x{d}), x{d} = full(x{d});   end

    X    = tensorprod(X,x{d}(:),DIM(d),1);
    DIM  = DIM - 1;
end

return

