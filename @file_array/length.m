function l = length(x)
% Overloaded length function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: length.m 174 2005-05-24 11:03:32Z john $


l = max(size(x));

