function res = fieldnames(this, varargin)
% Returns names of the fields in .other
% FORMAT res = fieldnames(this)
%
% An overloaded function...
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fieldnames.m 5025 2012-10-31 14:44:13Z vladimir $

res = fieldnames(this.other, varargin{:});