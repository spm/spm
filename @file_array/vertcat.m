function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% $Id$

o = cat(1,varargin{:});
return;
