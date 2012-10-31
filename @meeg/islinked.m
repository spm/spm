function out = islinked(this)
% True if the object is linked to data file
% FORMAT out = islinked(this)
% _________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: islinked.m 5025 2012-10-31 14:44:13Z vladimir $

this = check(this);
out  = isa(this.data, 'file_array');
