function spm_latex(c)
% Convert a job configuration structure into a series of LaTeX documents
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$

if nargin==0, c = spm_config; end;

fp = fopen('spm_manual.tex','w');
fprintf(fp,'\\documentclass[a4paper,titlepage]{book}\n');
fprintf(fp,'\\usepackage{epsfig,amsmath,pifont,moreverb,minitoc}\n');
fprintf(fp,'\\pagestyle{headings}\n\\bibliographystyle{authordate1}\n\n');
fprintf(fp,'\\hoffset=15mm\n\\voffset=-5mm\n');
fprintf(fp,'\\oddsidemargin=0mm\n\\evensidemargin=0mm\n\\topmargin=0mm\n');
fprintf(fp,'\\headheight=12pt\n\\headsep=10mm\n\\textheight=240mm\n\\textwidth=148mm\n');
fprintf(fp,'\\marginparsep=5mm\n\\marginparwidth=21mm\n\\footskip=10mm\n\n');

fprintf(fp,'\\title{\\huge{SPM5b (beta version) Manual}}\n');
fprintf(fp,'\\author{The FIL Methods Group}\n');
fprintf(fp,'\\begin{document}\n');
fprintf(fp,'\\maketitle\n');
fprintf(fp,'\\dominitoc\\tableofcontents\n\n');
% write_help(c,fp);
for i=1:numel(c.values),
   if isfield(c.values{i},'tag'),
       part(c.values{i},fp);
   end;
end;
fprintf(fp,'\\parskip=0mm\n\\bibliography{refs}\n\\end{document}\n\n');
fclose(fp);
return;

function part(c,fp)
if isstruct(c) && isfield(c,'tag'),
    fprintf(fp,'\\part{%s}\n',texify(c.name));
    % write_help(c,fp);
    if isfield(c,'values'),
        for i=1:numel(c.values),
            if isfield(c.values{i},'tag'),
                fprintf(fp,'\\include{%s}\n',c.values{i}.tag);
                chapter(c.values{i});
            end;
        end;
    end;
end;
return;

function chapter(c)
fp = fopen([c.tag '.tex'],'w');
fprintf(fp,'\\chapter{%s}\n\\minitoc\n\n',texify(c.name));
write_help(c,fp);

switch c.type,
case {'branch'},
    for i=1:numel(c.val),
        section(c.val{i},fp);
    end;
case {'repeat','choice'},
    for i=1:numel(c.values),
        section(c.values{i},fp);
    end;
end;
fclose(fp);
return;

function section(c,fp,lev)
if nargin<3, lev = 1; end;
sec = {'section','subsection','subsubsection','paragraph','subparagraph'};
if lev<=5,
    fprintf(fp,'\n\\%s{%s}\n',sec{lev},texify(c.name));
    write_help(c,fp);
    switch c.type,
    case {'branch'},
        for i=1:numel(c.val),
            section(c.val{i},fp,lev+1);
        end;
    case {'repeat','choice'},
        for i=1:numel(c.values),
            section(c.values{i},fp,lev+1);
        end;
    end;
else
    warning(['Too many nested levels... ' c.name]);
end;
return;

function write_help(hlp,fp)
if isstruct(hlp),
    if isfield(hlp,'help'),
       hlp = hlp.help;
    else
        return;
    end;
end;
if iscell(hlp),
    for i=1:numel(hlp),
        write_help(hlp{i},fp);
    end;
    return;
end;
str = texify(hlp);
fprintf(fp,'%s\n\n',str);
return;

function str1 = texify(str0)
[st,en,tok]=regexp(str0,'/\*([^(/\*)]*)\*/','start','end','tokens');
str1 = [];
st1  = [1 en+1];
en1  = [st-1 numel(str0)];
for i=1:numel(st1),
    str  = str0(st1(i):en1(i));
    str  = strrep(str,'$','\$');
    str  = strrep(str,'&','\&');
    str  = strrep(str,'_','\_');
   %str  = strrep(str,'\','$\\$');
    str  = strrep(str,'|','$|$');
    str  = strrep(str,'>','$>$');
    str  = strrep(str,'<','$<$');
    str1 = [str1 str];
    if i<numel(st1), str1 = [str1 tok{i}{1}]; end;
end;
return;
