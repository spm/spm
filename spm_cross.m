function [Y] = spm_cross(X,x)
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
% $Id: spm_cross.m 6812 2016-06-18 11:16:21Z karl $

% inner product using bsxfun
%--------------------------------------------------------------------------
A = reshape(full(X),[size(X) ones(1,ndims(x))]);
B = reshape(full(x),[ones(1,ndims(X)) size(x)]);
Y = squeeze(bsxfun(@times,A,B));
