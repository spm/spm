function [strout,R2] = spm_str_manip(strin,options)
% miscellaneous string manipulation options
% FORMAT string_out = spm_str_manip(string_in,options)
% string_in  - input string, string matrix, or cell array of strings
% string_out - output sring, string matrix, or cell array of strings
% options    - a string of options flags
%_______________________________________________________________________
% Each of the options is performed from left to right.
% The options are:
% 	'r'              - remove trailing suffix
% 	's'              - remove trailing suffix -
%			   only if it is either '.img', '.hdr', '.mat' or '.mnc'
% 	'e'              - remove everything except the suffix
% 	'h'              - remove trailing pathname component
% 	't'              - remove leading pathname component
% 	['f' num2str(n)] - remove all except first n characters
% 	['l' num2str(n)] - remove all except last n characters
%	['k' num2str(n)] - produce a string of at most n characters long.
%			   If the input string is longer than n, then
%			   it is prefixed with '..' and the last n-2 characters
%	['a' num2str(n)] - similar to above - except the leading directories
%			   are replaced by './'.
%			   eg. spm_str_manip('/dir1/dir2/file.img','a16') would
%			   produce '../dir2/file.img'.
%			   are returned.
%       'v'              - delete non valid filename characters
%                          Valid are './abc...xyzAbc...XYZ0123456789_-'
% 	'c'		 - remove leading components common to all strings
%                          returns leading component as a second output argument
%	'd'		 - deblank - this is always done!
%
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin<2, options=''; end
if nargin<1, strout=[]; R2=''; return, end

strout = cellstr(strin);
R2     = '';

while (~isempty(options))
	o = 2;

	% Read in optional numeric argument c
	%---------------------------------------------------------------
	opt = options(2:length(options));
	c = 0;
	if (~isempty(opt)), b = opt(1); else, b = '-'; end
	while(b >= '0'+0 & b <= '9'+0)
		c = c * 10 + opt(1)-'0';
		opt = opt(2:length(opt));
		if (isempty(opt))
			b = '-';
		else
			b = opt(1);
		end
		o = o + 1;
	end


	for i=1:prod(size(strout))
		str = deblank(strout{i});
		switch options(1)
		
		case 'r'	% Remove a trailing suffix of the form `.xxx',
				% leaving the basename.
			d1 = max([find(str == '/') 0]);
			d2 = max([find(str == '.') 0]);
			if (d2>d1), str = str(1:(d2-1)); end

		case 'e'	% Remove all but the suffix.
			d1 = max([find(str == '/') 0]);
			d2 = max([find(str == '.') 0]);
			if (d2>d1)
				str = str((d2+1):length(str));
			else
				str = '';
			end

		case 'h'	% Remove a trailing pathname component,
				% leaving the head.
			d1 = max([find(str == '/') 0]);
			if (d1>0)
				str = str(1:(d1-1));
			end

		case 't'	% Remove all leading  pathname  components,
				% leaving the tail.
			d1 = max([find(str == '/') 0]);
			if (d1>0)
				str = str((d1+1):length(str));
			end

		case 'f'	% First few characters
			str = str(1:min([length(str) c]));

		case 'l'	% Last few characters
			l = length(str);
			str = str(l-min([length(str) c])+1:l);

		case 'k'	% Last few characters
			l = length(str);
			if (l>c)
				str = ['..' str(l-c+2:l)];
			end

		case 'a'	% Last few characters
			m1   = find(str == '/');
			l    = length(str);
			if (c < l)
				m2   = find(l-m1+1+2 <= c);
				if ~isempty(m2)
					str = ['.' str(m1(min(m2)):l)];
				else
					str = ['.' str(max(m1):l)];
				end
			end

		case 's'	% Strip off '.img', '.hdr' or '.mat' suffixes
			l = length(str);
			if (l > 4)
				if (strcmp(str((l-3):l),'.img') | ...
				    strcmp(str((l-3):l),'.hdr') | ...
				    strcmp(str((l-3):l),'.mnc') | ...
				    strcmp(str((l-3):l),'.mat'))
					str = spm_str_manip(str, 'r');
				end
			end

		case 'v'
			tmp = find(...
				( str>='a' & str<='z' ) | ...
				( str>='A' & str<='Z' ) | ...
				( str>='0' & str<='9' ) | ...
				  str=='-' | str=='_' | ...
				  str=='.' | str=='/' );
			str = str(tmp);

		end
		strout{i} = str;
	end
		
	if options(1)=='c'	% Remove common path components
		if length(strout)>1
			tmp    = size(strout);
			strout = char(strout(:));
			msk    = diff(strout+0)~=0;
			d1     = min(find(sum(msk)));
			d1     = max([find(strout(1,1:d1) == '/') 0]);
			R2     = strout(1,1:d1);
			strout = reshape(cellstr(strout(:,d1+1:end)),tmp);
		end

	elseif options(1)=='d'	% Deblanking is always done by cellstr above
		% strout = deblank(strout);
	end

	options = options(o:length(options));	

end

if ischar(strin), strout=char(strout); end
