function out = cfg_justify(varargin)
% CFG_JUSTIFY Justifies a text string
%    OUT = CFG_JUSTIFY(N,TXT) justifies text string TXT to
%    the length specified by N.
%
%    OUT = CFG_JUSTIFY(OBJ,TXT), where OBJ is a handle to a 'listbox' style
%    uicontrol, justifies text string TXT to the width of the OBJ in
%    characters - 1.
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
%    out = cfg_justify(40,{['Statistical Parametric ',...
%    'Mapping refers to the construction and ',...
%    'assessment of spatially extended ',...
%    'statistical process used to test hypotheses ',...
%    'about [neuro]imaging data from SPECT/PET & ',...
%    'fMRI. These ideas have been instantiated ',...
%    'in software that is called SPM']});
%    strvcat(out{:})
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: cfg_justify.m 3944 2010-06-23 08:53:40Z volkmar $

out = {};

if nargin < 2
    cfg_message('matlabbatch:usage','Incorrect usage of cfg_justify.')
end
n = varargin{1};
if ishandle(n)
    % estimate extent of a space char and scrollbar width
    TempObj=copyobj(n,get(n,'Parent'));
    set(TempObj,'Visible','off','Max',100);
    spext = cfg_maxextent(TempObj,{repmat(' ',1,100)})/100;
    % try to work out slider size
    pos = get(TempObj,'Position');
    oldun = get(TempObj,'units');
    set(TempObj,'units','points');
    ppos = get(TempObj,'Position');
    set(TempObj,'units',oldun);
    sc = pos(3)/ppos(3);
    % assume slider width of 15 points
    swidth=15*sc;
else
    % dummy constants
    spext = 1;
    swidth = 0;
    TempObj = n;
end
for i=2:nargin,
    if iscell(varargin{i}),
        for j=1:numel(varargin{i}),
            para = justify_paragraph(TempObj,spext,swidth,varargin{i}{j});
            out  = [out(:);para(:)]';
        end
    else
        para = justify_paragraph(TempObj,spext,swidth,varargin{i});
        out  = [out(:);para(:)]';
    end
end
if ishandle(TempObj)
    delete(TempObj);
end
function out = justify_paragraph(n,spext,swidth,txt)
if numel(txt)>1 && txt(1)=='%',
    txt = txt(2:end);
end;
%txt = regexprep(txt,'/\*([^(/\*)]*)\*/','');
st1  = strfind(txt,'/*');
en1  = strfind(txt,'*/');
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
if isempty(off),
    out{1} = txt;
else
    off  = off(1);
    para = justify_para(n,off,spext,swidth,txt(off:end));
    out = cell(numel(para),1);
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

function out = justify_para(n,off,spext,swidth,varargin)
% Collect varargs into a single string
str = varargin{1};
for i=2:length(varargin),
    str = [str ' ' varargin{i}];
end;
if isempty(str), out = {''}; return; end;
if ishandle(n)
    % new size: at least 20 spaces wide, max widget width less offset
    % space and scrollbar width
    pos = get(n,'position');
    opos = pos;
    pos(3) = max(pos(3)-off*spext-swidth,20*spext);
    set(n,'position',pos);
    % wrap text
    out = textwrap(n,{str});
    cext = cfg_maxextent(n,out);
    % fill with spaces to produce (roughly) block output
    for k = 1:numel(out)-1
        out{k} = justify_line(out{k}, pos(3), cext(k), spext);
    end;
    % reset position
    set(n,'Position',opos);
else
    cols = max(n-off,20);
    out = textwrap({str},cols);
    for k = 1:numel(out)-1
        out{k} = justify_line(out{k}, cols, length(out{k}), 1);
    end;
end;

function out = justify_line(str, width, cext, spext)
ind = strfind(str,' ');
if isempty(ind)
    out = str;
else
    % #spaces to insert
    nsp = floor((width-cext)/spext);
    % #spaces per existing space
    ins(1:numel(ind)) = floor(nsp/numel(ind));
    ins(1:mod(nsp,numel(ind))) = floor(nsp/numel(ind))+1;
    % insert spaces beginning at the end of the string
    for k = numel(ind):-1:1
        str = [str(1:ind(k)) repmat(' ',1,ins(k)) str(ind(k)+1:end)];
    end;
    out = str;
end; 
