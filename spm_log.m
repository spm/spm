function LogFile = spm_log(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)
% SPM logging function. Writes string arguments out to a log file.
% FORMAT LogFile = spm_log(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)
% P{1-10} - Matrices to be written to the log file.
% LogFile - The name of the log file.
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
global LOGFILE; LogFile=LOGFILE;
if isempty(LogFile), return, end

if nargin==0, return, end

%-Open LogFile
%----------------------------------------------------------------------------
[fid,message] = fopen(LogFile,'a');
if fid==-1
	fprintf('%cspm_log error : tNo logging!\n\t%s\n\n',7,message);
	return
end


%-Write log
%----------------------------------------------------------------------------
for arg = 1:nargin
	tmp = eval(['P',int2str(arg)]);
	if isstr(tmp)
		Str = tmp;
	else
		%-Build string matrix representation of numeric matrix
		fmt = '%8.6g ';
		% if all(floor(tmp(:))==ceil(tmp(:)))
		%	fmt = ['%',int2str(ceil(log10(max([1,abs(tmp(:))])))+1),'d ']; end
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
