function a=subsasgn(a,b,c)
% Overloaded subsasgn function for spm_file_array objects.
% _______________________________________________________________________
% @(#)subsasgn.m	1.1 John Ashburner 04/06/28
error('Sorry, but the elements of spm_file_array objects are read-only.');
