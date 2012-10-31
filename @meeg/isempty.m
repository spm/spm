function out = isempty(this)
% True if the object is empty 
% FORMAT out = isempty(this)
% _________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: isempty.m 5025 2012-10-31 14:44:13Z vladimir $

out = all(size(this)==0);
