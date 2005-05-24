function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: vertcat.m 174 2005-05-24 11:03:32Z john $


o = cat(1,varargin{:});
return;
