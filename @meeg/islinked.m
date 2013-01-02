function out = islinked(this)
% True if the object is linked to data file
% FORMAT out = islinked(this)
% _________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: islinked.m 5167 2013-01-02 15:24:52Z vladimir $

out  = isa(this.data, 'file_array');
