function varargout = spm_unlink(varargin)
% Silently delete files on disk - a compiled routine
% FORMAT spm_unlink file1 file2 file3 file4...
%     OR spm_unlink('file1','file2','file3','file4',...)
%
%__________________________________________________________________________
% Copyright (C) 1996-2020 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_unlink.m 7833 2020-04-17 10:43:06Z guillaume $


%-This is merely the help file for the compiled routine
%error('spm_unlink.c not compiled - see Makefile')

rs  = recycle('off');
crs = onCleanup(@() recycle(rs));
ws  = warning('off');
wsr = onCleanup(@() warning(ws));
for i=1:numel(varargin)
    delete(varargin{i});
end
