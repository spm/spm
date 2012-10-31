function this = unlink(this)
% Unlinks the object from the data file 
% FORMAT this = unlink(this)   
% _________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: unlink.m 5025 2012-10-31 14:44:13Z vladimir $

this.data = [];
this      = check(this);