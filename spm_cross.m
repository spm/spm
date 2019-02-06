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
% $Id: spm_cross.m 7527 2019-02-06 19:12:56Z karl $

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
A = reshape(full(X),[size(X) ones(1,ndims(x))]);
B = reshape(full(x),[ones(1,ndims(X)) size(x)]);
Y = squeeze(bsxfun(@times,A,B));

% and handle remaining arguments
%--------------------------------------------------------------------------
for i = 1:numel(varargin)
    Y = spm_cross(Y,varargin{i});
end
