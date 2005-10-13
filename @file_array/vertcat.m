function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: vertcat.m 253 2005-10-13 15:31:34Z guillaume $


o = cat(1,varargin{:});
return;
