function l = length(x)
% Overloaded length function for file_array objects
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: length.m 7147 2017-08-03 14:07:01Z spm $


l = max(size(x));
