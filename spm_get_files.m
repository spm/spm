function varargout = spm_get_files(varargin)
% Extended gateway to spm_list_files: Returns full pathnames & canonicalises dir
% FORMAT Q spm_get_files(dir,fil)
% dir      - directory within which list files (a string)
%            parameter is canonicalised by spm_get('CPath',dir)
%            [defaults to '.', the current directory]
% fil      - file filter string: e.g. 'sn*.img' [default '*']
% P        - string matrix of full pathnames to files in directory
% dir      - full pathname of directory listed (after canonicalisation)
%_______________________________________________________________________
%
% spm_results functionality has been absorbed into spm_get.
%
%                           ----------------
%
% Based on spm_list_files.m
% Used in m-files for batch mode to help entering the file names.
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_get_files is grandfathered, use spm_get(''Files'',dir,fil) instead')


%-Pass on to spm_get & return necessary arguments
%-----------------------------------------------------------------------
varargout = {spm_get('Files',varargin{:})};
