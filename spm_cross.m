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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cross.m 7838 2020-04-23 17:40:45Z karl $

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
Y   = reshape(Y,[siz(siz > 1) 1]);

% and handle remaining arguments
%--------------------------------------------------------------------------
for i = 1:numel(varargin)
    Y = spm_cross(Y,varargin{i});
end
