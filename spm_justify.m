function out = spm_justify(n,varargin)
% Justify text
% FORMAT out = justify(n,txt)
% out - a cell array of lines of text
% n   - line length
% txt - a text string or a cell array of text strings
%
% If txt is a cell array, then each element is treated
% as a paragraph and justified, otherwise the string is
% treated as a paragraph and is justified.
% Non a-z or A-Z characters at the start of a paragraph
% are used to define any indentation required (such as
% for enumeration, bullets etc.  If less than one line
% of text is returned, then no formatting is done.
%
% Example usage:
%
%    out = spm_justify(40,{['Statistical Parametric ',...
%    'Mapping refers to the construction and ',...
%    'assessment of spatially extended ',...
%    'statistical process used to test hypotheses ',...
%    'about [neuro]imaging data from SPECT/PET & ',...
%    'fMRI. These ideas have been instantiated ',...
%    'in software that is called SPM']});
%    strvcat(out{:})
%
%------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$

out = {};
for i=1:nargin-1,
    if iscell(varargin{i}),
        for j=1:numel(varargin{i}),
            para = justify_paragraph(n,varargin{i}{j});
            out  = {out{:},para{:}};
        end;
    else
        para = justify_paragraph(n,varargin{i});
        out = {out{:},para{:}};
    end;
end;

function out = justify_paragraph(n,txt)
txt = regexprep(txt,'/\*([^(/\*)]*)\*/','');
off = find((txt'>='a' & txt'<='z') | (txt'>='A' & txt'<='Z'));
off = off(off<n);
if isempty(off),
    out{1} = txt;
else
    off  = off(1);
    para = justify_para(n-off+1,txt(off:end));
    if numel(para)>1,
        out{1} = [txt(1:(off-1)) para{1}];
        for j=2:numel(para),
            out{j} = [repmat(' ',1,off-1) para{j}];
        end;
    else
        out{1} = txt;
    end;
end;
return;

function out = justify_para(n,varargin)
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
