function l = length(x)
% Overloaded length function for spm_file_array objects
%
% _______________________________________________________________________
% @(#)length.m	1.1 John Ashburner 04/06/28

l = max(size(x));

