function LOGFILE = spm_log(varargin)
% SPM logging function. Writes string arguments out to a log file.
% FORMAT LOGFILE = spm_log(varargin)
% varargin - Matrices to be written to the log file.
% LOGFILE  - The name of the log file.
%__________________________________________________________________________
%
% spm_log implements logging for the SPM package. The log file is
% specified in the global LOGFILE. If this is non-empty, then spm_log
% writes the passed string matrices to the log file.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%


%-Find out LogFile name, return if not logging.
%============================================================================
if any(strcmp(who('global'),'LOGFILE')), global LOGFILE; else, LOGFILE=''; end
if nargin==0 | isempty(LOGFILE), return, end

%-Open LogFile
%----------------------------------------------------------------------------
[fid,msg] = fopen(LOGFILE,'a');
if fid==-1
    if strcmp(msg,'Sorry. No help in figuring out the problem . . .')
	warning(sprintf('spm_log error:  No write permission for %s',...
		spm_str_manip(tmp,'p')))
    else    
	warning(['spm_log error: ',msg])
    end    
    return
end


%-Write log
%----------------------------------------------------------------------------
for arg = 1:nargin
	tmp = varargin{arg};
	if isstr(tmp)
		Str = tmp;
	elseif isempty(tmp)
	        Str = '[]';
	else
		%-Build string matrix representation of numeric matrix
		Str = '';
		for r = 1:size(tmp,1)
			Str = strvcat(Str,sprintf('%8.6g ',tmp(r,:))); end
	end
	for str = Str', fprintf(fid,[str','\n']); end
end
fprintf(fid,'\n');

%-Close log
%----------------------------------------------------------------------------
status = fclose(fid);
if fid==-1, warning('spm_log : error closing file'), end
