function res = transformtype(this, varargin)
% Method for getting/setting type of transform
% FORMAT res = transformtype(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: transformtype.m 1270 2008-03-28 14:35:16Z stefan $

res = getset(this, 'transform', 'ID', 1, varargin{:});
