function strout = spm_str_manip(strin, options)
% miscellaneous string manipulation options
% FORMAT string_out = spm_str_manip(string_in,options)
% string_in  - input string (or string matrix)
% string_out - output sring (or string matrix)
% options    - a string of options flags
%____________________________________________________________________________
% Each of the options is performed from left to right.
% The options are:
% 	'r'              - remove trailing suffix
% 	's'              - remove trailing suffix -
%			   only if it is either '.img' or '.hdr'
% 	'e'              - remove everything except the suffix
% 	'h'              - remove trailing pathname component
% 	't'              - remove leading pathname component
% 	['f' num2str(n)] - remove all except first n characters
% 	['l' num2str(n)] - remove all except last n characters
%	['k' num2str(n)] - produce a string of at most n characters long.
%			   If the input string is longer than n, then
%			   it is prefixed with '..' and the last n-2 characters
%			   are returned.
% 	'c'		 - remove leading components which are common to
%			   all rows of string-matrix
%	'd'		 - deblank
%
%____________________________________________________________________________
% %W% John Ashburner %E%

strout = strin;

while (~isempty(options))
	m = size(strout,1);
	n = size(strout,2);
	o = 2;


	% Read in optional numeric argument c
	opt = options(2:length(options));
	c = 0;
	b = '-';
	if (~isempty(opt))b = opt(1);end

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


	for i=1:m
		opt = options;
		stri = deblank(strout(i,:));
		if (opt(1)=='r')
			% Remove a trailing suffix of the form `.xxx',
			% leaving the basename.
			d1 = max([find(stri == '/') 0]);
			d2 = max([find(stri == '.') 0]);
			if (d2>d1), stri = stri(1:(d2-1)); end

		elseif (opt(1)=='e')
			% Remove all but the suffix.
			d1 = max([find(stri == '/') 0]);
			d2 = max([find(stri == '.') 0]);
			if (d2>d1)
				stri = stri((d2+1):length(stri));
			else
				stri = '';
			end

		elseif (opt(1)=='h')
			% Remove a trailing pathname component,
			% leaving the head.
			d1 = max([find(stri == '/') 0]);
			if (d1>0)
				stri = stri(1:(d1-1));
			end

		elseif (opt(1)=='t')
			% Remove all leading  pathname  components,
			% leaving the tail.
			d1 = max([find(stri == '/') 0]);
			if (d1>0)
				stri = stri((d1+1):length(stri));
			end

		elseif (opt(1)=='f')
			% First few characters
			stri = stri(1:min([length(stri) c]));

		elseif (opt(1)=='l')
			% Last few characters
			l = length(stri);
			stri = stri(l-min([length(stri) c])+1:l);

		elseif (opt(1)=='k')
			% Last few characters
			l = length(stri);
			if (l>c)
				stri = ['..' stri(l-c+2:l)];
			end

		elseif (opt(1)=='s')
			% Strip off '.img' or '.hdr' suffixes only
			l = length(stri);
			if (l > 4)
				if (strcmp(stri((l-3):l),'.img') | strcmp(stri((l-3):l),'.hdr'))
					stri = spm_str_manip(stri, 'r');
				end
			end
		end
		l = length(stri);
		if (l<n)
			strout(i,:) = [stri zeros(1,(n-l))+' '];
		else
			strout(i,:) = stri;
		end

	end
		
	if (options(1) == 'c')
		% Remove common path components..
		if (m>1)
			msk = diff(strout+0) ~=0;
			d1 = min(find(sum(msk)));
			d1 = max([find(strout(1,1:d1) == '/') 0])+1;
			strout = [strout(:,d1:n) zeros(m,(d1-1))+' '];
		end


	elseif (options(1) == 'd')
		% Deblank
		if (m<=1)
			strout = deblank(strout);
		else
			msk = fliplr(sum(strout~=' '));
			d1  = min(find(msk));
			strout = strout(:,1:(n-d1+1));
		end
	end

	options = options(2:length(options));	

end

