function SPMdir = spm_get_path
% returns the SPM path or SPM directory
% FORMAT SPMdir = spm_get_path
%___________________________________________________________________________
%
% spm_get_path returns the SPM95 directory
%
% This directory is found as the first on MatLab's path containing
% the SPM95 startup m-file "spm".
%
% This function is called by spm_ui which sets the global variable SWD
% (SPM Working Directory) to this directory's pathname
%
% If this function generates an error, then the SPM directory is not on
% MatLabs path.
%
% To append the SPM directory (e.g. /usr/brain/spm) to MATLABPATH
% add the following to your .cshrc and then login again
%
% setenv MATLABPATH "/usr/brain/spm:$MATLABPATH"
%
%__________________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%


% Version History
% - Karl Friston  - V1 - 04/94 - Searched for file .spm94 in every
%                                directory on the path
% - Andrew Holmes - V2 - 02/95 - Improved speed by using built in
%                                "which" command


%-Computation
%=============================================================================

Mfile = 'spm';

str = which(Mfile);

if ~isstr(str)
	error('SPM directory not in MATLABPATH:  type help spm_get_path'); end

SPMdir = strrep(str,['/',Mfile,'.m'],'');

if isempty(SPMdir), SPMdir='/'; end
