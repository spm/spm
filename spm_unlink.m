function spm_unlink(varargin)
% Silently delete files on disk - a compiled routine
% FORMAT spm_unlink('file1','file2','file3','file4',...)
%
% Remove the specified file(s) using a system call to unlink().
%__________________________________________________________________________
% Copyright (C) 1996-2020 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_unlink.m 7840 2020-04-26 23:11:25Z spm $


%-This is merely the help file for the compiled routine
%error('spm_unlink.c not compiled - see Makefile')

rs  = recycle('off');
crs = onCleanup(@() recycle(rs));

ws  = warning('off');
cws = onCleanup(@() warning(ws));

for i=1:numel(varargin)
    delete(varargin{i});
end
