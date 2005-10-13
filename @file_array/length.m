function l = length(x)
% Overloaded length function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: length.m 253 2005-10-13 15:31:34Z guillaume $


l = max(size(x));

