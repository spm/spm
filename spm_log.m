function LOGFILE = spm_log(varargin)
% SPM logging function. Writes string arguments out to a log file.
% FORMAT LOGFILE = spm_log(varargin)
% varargin - Matrices to be written to the log file.
% LOGFILE  - The name of the log file.
%_______________________________________________________________________
%
% spm_log implements logging for the SPM package.
%
% The log file is specified in the global LOGFILE.
%
% If this is non-empty, then spm_log writes the passed string matrices 
% to the log file.
%
% Doesn't handle multi-dimensional (>2) matrices very gracefully!
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%


%-Find out LogFile name, return if not logging.
%=======================================================================
LOGFILE = spm('GetGlobal','LOGFILE');
if nargin==0 | isempty(LOGFILE), return, end

%-Open LogFile
%-----------------------------------------------------------------------
[fid,msg] = fopen(LOGFILE,'a');
if fid==-1
	if strcmp(msg,'Sorry. No help in figuring out the problem . . .')
		msg = 'No write permission';
	end    
	spm('alert!',{	'Problems logging input to LOGFILE:',' ',...
		['      ',spm_str_manip(LOGFILE,'p')],' ',...
		['-> ',msg],' ',...
		'Resetting LOGFILE'},mfilename)
	global LOGFILE, LOGFILE = '';
	return
end


%-Write log
%-----------------------------------------------------------------------
for arg = 1:nargin
	if isstr(varargin{arg})
		Str = cellstr(varargin{arg});
	elseif isempty(varargin{arg})
	        Str = {'[]'};
	elseif iscellstr(varargin{arg})
		Str = varargin{arg};
	else
		%-Build string matrix representation of numeric matrix
		tmp = varargin{arg};
		for r = 1:size(tmp,1), Str{r} = sprintf('%8.6g ',tmp(r,:)); end
	end
	fprintf(fid,'%s\n',Str{:});
end
fprintf(fid,'\n');

%-Close log
%-----------------------------------------------------------------------
status = fclose(fid);
if fid==-1
	spm('alert!',{	'Error closing LOGFILE:',' ',...
		['      ',spm_str_manip(LOGFILE,'p')],' ',...
		'Resetting LOGFILE'},mfilename)
	global LOGFILE, LOGFILE = '';
end
