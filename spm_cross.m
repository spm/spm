function [Y] = spm_cross(X,x,varargin)
% Multidimensional cross (outer) product
% FORMAT [Y] = spm_cross(X,x)
% FORMAT [Y] = spm_cross(X)
%
% X  - numeric array
% x  - numeric array
%
% Y  - outer product
%
% See also: spm_dot
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging


% handle single inputs
%--------------------------------------------------------------------------
if nargin < 2
    if isnumeric(X)
        Y = X;
    else
        Y = spm_cross(X{:});
    end
    return
end

% handle cell arrays
%--------------------------------------------------------------------------
if iscell(X), X = spm_cross(X{:}); end
if iscell(x), x = spm_cross(x{:}); end

% outer product of first pair of arguments (using bsxfun)
%--------------------------------------------------------------------------
A   = reshape(full(X),[size(X) ones(1,ndims(x))]);
B   = reshape(full(x),[ones(1,ndims(X)) size(x)]);
Y   = bsxfun(@times,A,B);
siz = size(Y);
siz = siz(siz > 1);
if numel(siz)
    Y = reshape(Y,[siz,1]);
end

% and handle remaining arguments
%--------------------------------------------------------------------------
for i = 1:numel(varargin)
    Y = spm_cross(Y,varargin{i});
end
