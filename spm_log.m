function LogFile = spm_log(varargin)
% SPM logging function. Writes string arguments out to a log file.
% FORMAT LogFile = spm_log(varargin)
% varargin - Matrices to be written to the log file.
% LogFile  - The name of the log file.
%__________________________________________________________________________
%
% spm_log implements logging for the SPM package. The log file is
% specified in the global LOGFILE. If this is non-empty, then spm_log
% writes the passed string matrices to the log file.
%
%__________________________________________________________________________
% @(#)spm_log.m	1.2 Andrew Holmes 97/08/08


%-Find out LogFile name, return if not logging.
%============================================================================
global LOGFILE; LogFile=LOGFILE;
if isempty(LogFile), return, end

if nargin==0, return, end

%-Open LogFile
%----------------------------------------------------------------------------
[fid,message] = fopen(LogFile,'a');
if fid==-1
    if strcmp(message,'Sorry. No help in figuring out the problem . . .')
	tmp = LogFile;
	if (tmp(1)~='/'); tmp = [pwd '/' tmp]; end
	fprintf('%cspm_log error:  No write permission for ''%s''\n',...
		7,tmp);
    else    
	fprintf('%cspm_log error:  No logging!\n\t%s\n\n',7,message);
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
		fmt = '%8.6g ';
		%if all(floor(tmp(:))==ceil(tmp(:)))
		%fmt=['%',int2str(ceil(log10(max([1,abs(tmp(:))])))+1),'d '];end
		Str = sprintf(fmt,tmp(1,:));
		for r = 2:size(tmp,1)
			Str = str2mat(Str,sprintf(fmt,tmp(r,:))); end
	end
	for str = Str', fprintf(fid,[str','\n']); end
end
fprintf(fid,'\n');

%-Close log
%----------------------------------------------------------------------------
status = fclose(fid);
if fid==-1, fprintf('%cspm_log : error closing file, fid = %d\n',7,fid); end
