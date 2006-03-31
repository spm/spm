function out = spm_justify(varargin)
% SPM_JUSTIFY Justifies a text string
%    OUT = SPM_JUSTIFY(N,TXT) justifies text string TXT to
%    the length specified by N.
%
%    OUT = SPM_JUSTIFY(OBJ,TXT) justifies text string to
%    the width of the OBJ in characters - 1.
%
%    If TXT is a cell array, then each element is treated
%    as a paragraph and justified, otherwise the string is
%    treated as a paragraph and is justified.
%    Non a-z or A-Z characters at the start of a paragraph
%    are used to define any indentation required (such as
%    for enumeration, bullets etc.  If less than one line
%    of text is returned, then no formatting is done.
%
%    Example:
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
% $Id: spm_justify.m 491 2006-03-31 17:36:45Z john $

out = {};

% DRG added code to estimate the actual size of the listbox and how many
% characters will fit per line. 
if nargin < 2
    error('Incorrect usage of spm_justify.')
end
n = varargin{1};
if ishandle(n)
    % n is presumably the handle of some object that can deal with text.
    % use try catch statements so that if try fails the
    % UI is returned to its original state.
    h = n;
    oldStr = get(h,'String');
    oldUnits = get(h,'Units');
    try
        set(h,'Units','Characters');
        % Matlab is very inaccurate at estimating the width of single
        % characters, so we'll use 100.
        set(h,'string',repmat(' ',1,100));
        ext = get(h,'extent');
        % calculate width of 1 character.
        ext = ext(3)/100;
        pos = get(h,'position');
        pos = pos(3);
        % hack because Matlab cannot tell us how wide the scrollbar is.
        pos = pos - 5;
        wGuess = floor(pos/ext);
        str = repmat(' ',1,wGuess);
        set(h,'string',str);
        ext = get(h,'extent');
        ext = ext(3);
        if ext > pos
            while ext > pos
                if numel(str) == 1
                    error('Too few characters for object handle provided.');
                end
                str = str(1:end-1);
                set(h,'String',str);
                ext = get(h,'extent');
                ext = ext(3);
            end
            n = numel(str);
        else
            while ext < pos
                str = [str(1), str];
                if numel(str) > 200
                    error('Too many characters for object handle provided.');
                end
                set(h,'String',str);
                ext = get(h,'extent');
                ext = ext(3);
            end
            n = numel(str)-2;
        end
        set(h,'String',oldStr);
        set(h,'Units',oldUnits);
    catch
        set(h,'String',oldStr);
        set(h,'Units',oldUnits);
        error('Problem estimating characters for object handle provided');
    end
end
for i=2:nargin,
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
if numel(txt)>1 && txt(1)=='%',
    txt = txt(2:end);
end;
%txt = regexprep(txt,'/\*([^(/\*)]*)\*/','');
st1  = findstr(txt,'/*');
en1  = findstr(txt,'*/');
st = [];
en = [];
for i=1:numel(st1),
    en1  = en1(en1>st1(i));
    if ~isempty(en1),
        st  = [st st1(i)];
        en  = [en en1(1)];
        en1 = en1(2:end);
    end;
end;

str = [];
pen = 1;
for i=1:numel(st),
    str = [str txt(pen:st(i)-1)];
    pen = en(i)+2;
end;
str = [str txt(pen:numel(txt))];
txt = str;

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
    % l   = length(current);
    % l   = n-l;
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
        [unused,ind] = sort(rand(pad,1));
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
