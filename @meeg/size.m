function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 3210 2009-06-17 13:46:25Z vladimir $

res = size(this.data.y, varargin{:});
