function out = spm_justify(n,varargin)
% Justify text
% FORMAT out = justify(n,varargin)
% out - a cell array of lines of text
% n   - line length
% varargin - various bits of a paragraph
%
% Example usage:
%
%    out = spm_justify(40,'Statistical Parametric',...
%    'Mapping refers to the construction and',...
%    'assessment of spatially extended',...
%    'statistical process used to test hypotheses',...
%    'about [neuro]imaging data from SPECT/PET &',...
%    'fMRI. These ideas have been instantiated',...
%    'in software that is called SPM');
%    strvcat(out{:})
%
%------------------------------------------------------------------------
% John Ashburner $Id$


% Collect varargs into a single string
str = varargin{1};
for i=2:length(varargin),
	str = [str ' ' varargin{i}];
end;

if isempty(str), out = {''}; return; end;

% Remove repeats
sp  = find(str==' ');
rep = sp(diff(sp)==1);
str(rep) = [];
if str(1)  ==' ', str(1)   = ''; end;
if str(end)==' ', str(end) = ''; end;

out = {};
while length(str)>n,

	% Break the string into lines
	sp  = find(str==' ');
	brk = sp(sp<=n);
	if isempty(brk),
		if isempty(sp),
			brk = length(str)+1;
        else
			brk = sp(1);
		end;
    else
		brk = brk(end);
	end;

	% Pad the line to n characters wide
	current = str(1:(brk-1));
	l   = length(current);
	l   = n-l;
	sp  = find(current==' ');
	if ~isempty(sp),

		% Break into words
		sp    = [sp length(current)+1];
		words = {current(1:(sp(1)-1))};
		for i=1:(length(sp)-1),
			words = {words{:}, current((sp(i)+1):(sp(i+1)-1))};
		end;

		% Figure out how much padding on average
		nsp = length(sp)-1;
		pad = (n-(length(current)-nsp))/nsp;

		% Pad all spaces by the same integer amount
		sp  = repmat(floor(pad),1,nsp);

		% Pad a random selection of spaces by one
		pad = round((pad-floor(pad))*nsp);
		r   = rand(pad);
		[r,ind] = sort(rand(pad,1));
		ind = ind(1:pad);
		sp(ind) = sp(ind)+1;

		% Re-construct line from individual words
		current = words{1};
		for i=2:length(words),
			current = [current repmat(' ',1,sp(i-1)) words{i}];
		end;
	end;

	out = {out{:},current};
	str = str((brk+1):end);
end;

out = {out{:},str};
