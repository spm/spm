function s = spm_existfile(filename)
% Check if a file exists on disk - a compiled routine
% FORMAT S = SPM_EXISTFILE(FILENAME)
% FILENAME - filename (can also be a relative or full pathname to a file)
% S        - logical scalar, true if the file exists and false otherwise
%_______________________________________________________________________
%
% This compiled routine is equivalent to:
% >> s = exist(filename,'file') == 2;
% ans was written for speed purposes. The differences in behaviour are:
%  * spm_existfile returns true for directory names
%  * spm_existfile does not look in MATLAB's search path
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_existfile.m 1801 2008-06-09 19:00:47Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_existfile.c not compiled - see Makefile')
