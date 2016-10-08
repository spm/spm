function [Y] = spm_cross(X,x,varargin)
% Multidimensional cross (outer) product
% FORMAT [Y] = spm_cross(X,x)
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
% $Id: spm_cross.m 6901 2016-10-08 13:21:41Z karl $

% handle single inputs
%--------------------------------------------------------------------------
if nargin < 2, x = 1; end

% handle multiple inputs
%--------------------------------------------------------------------------
for i = 1:numel(varargin)
    X  = spm_cross(X,varargin{i});
end


% handle cell arrays
%--------------------------------------------------------------------------
if iscell(X), X = spm_cross(X{:}); end
if iscell(x), x = spm_cross(x{:}); end

% inner product using bsxfun
%--------------------------------------------------------------------------
A = reshape(full(X),[size(X) ones(1,ndims(x))]);
B = reshape(full(x),[ones(1,ndims(X)) size(x)]);
Y = squeeze(bsxfun(@times,A,B));
