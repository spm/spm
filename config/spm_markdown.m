function spm_markdown(c)
% Convert a job configuration tree into a series of markdown documents
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if ~nargin, c = spm_cfg; end
if nargin && ischar(c), clean_latex_compile; return; end

%cd(fullfile(spm('Dir'),'man'));

fp = fopen(fullfile('docs','spm_manual.md'),'w');
sectioning(c,fp);
bibcstr = get_bib(fullfile(spm('dir'),'man','biblio'));

tbxlist = dir(fullfile(spm('dir'),'toolbox'));
for k = 1:numel(tbxlist)
    if tbxlist(k).isdir
        bibcstr=[bibcstr(:); get_bib(fullfile(spm('dir'),'toolbox', ...
                 tbxlist(k).name))];
    end
end
bibcstr = strcat(bibcstr,',');
bibstr  = strcat(bibcstr{:});
fclose(fp);


function sectioning(c,fp,level,dr)
if nargin<3, level = 1; end
if nargin<4, dr = ''; end
switch class(c)
    case {'cfg_branch','cfg_exbranch'}
        if isa(c.val,'function_handle'), c.val = feval(c.val); end
end
switch class(c)
    case {'cfg_exbranch'}
        fprintf(fp,'\n\n%s [**%s**](../%s%s/)  \n',sec(level),texify(c.name),dr,c.tag);
        chapter(c, [dr c.tag '.md']);

    case {'cfg_branch'}
        fprintf(fp,'\n\n%s **%s**  \n', sec(level), texify(c.name));
        write_help(c,fp,level);

    case {'cfg_repeat','cfg_choice','cfg_menu'}
        fprintf(fp,'\n\n%s **%s**  \n', sec(level), texify(c.name));
        write_help(c,fp,level);
        for i=1:numel(c.values)
            sectioning(c.values{i},fp, level+1, [dr c.tag '.']);
        end

end

%==========================================================================
function s = sec(n)
if nargin<1, n=1; end
s = [repmat('    ',[1 n-1]) '*'];


%==========================================================================
function sts = chapter(c,fname)
if nargin<2, fname = c.tag; end
fp = fopen(fullfile('docs',fname),'w');
if fp==-1, sts = false; return; end

fprintf(fp,'# %s  \n', texify(c.name));
write_help(c,fp);

switch class(c)
    case {'cfg_branch','cfg_exbranch'}
        if isa(c.val,'function_handle'), c.val = feval(c.val); end
        for i=1:numel(c.val)
            section(c.val{i},fp);
        end
    case {'cfg_repeat','cfg_choice'}
        for i=1:numel(c.values)
            section(c.values{i},fp);
        end
end
fclose(fp);
sts = true;

%==========================================================================
function section(c,fp,lev)
if nargin<3, lev = 1; end

if lev>6 || strcmp(c.name,'Generic')
    return
else
    switch class(c)
        case 'cfg_choice'
            t = '(choose an option)';
        case 'cfg_const'
            t = '';
        case 'cfg_exbranch'
            t = '';
        case 'cfg_branch'
            t = '';
        case 'cfg_repeat'
            t = '(create a list of items)';
        case 'cfg_menu'
            t = '(choose from the menu)';
        case 'cfg_entry'
            t = '(enter text)';
        case 'cfg_files'
            t = '(select files)';
        otherwise
            t = ['(' class(c) ')'];
    end

    fprintf(fp,'\n%s **%s** %s  \n',sec(lev),texify(c.name), t);
    write_help(c,fp, lev);
    switch class(c)
        case {'cfg_branch','cfg_exbranch'}
            if isa(c.val,'function_handle'), c.val = feval(c.val); end
            for i=1:numel(c.val)
                section(c.val{i},fp,lev+1);
            end
        case {'cfg_repeat','cfg_choice'}
            for i=1:numel(c.values)
                section(c.values{i},fp,lev+1);
            end
    end
end


%==========================================================================
function write_help(hlp,fp,level)
if nargin<3, level=1; end
if isa(hlp, 'cfg_item')
    if ~isempty(hlp.help)
       hlp = hlp.help;
    else
        return;
    end
end
if iscell(hlp)
    for i=1:numel(hlp)
        if ~(numel(hlp{i}>1) && hlp{i}(1)=='%')
            write_help(hlp{i},fp, level);
        end
    end
    return;
end
str    = texify(hlp);
indent = repmat('    ',[1 level-1]);
if ~isempty(str)
    fprintf(fp,'%s%s  \n', indent, str);
else
    if true %indent>=1
        fprintf(fp,'%s.  \n', indent);
    else
        fprintf(fp,'  \n');
    end
end


%==========================================================================
function str = texify(str0)
st1  = strfind(str0,'/*');
en1  = strfind(str0,'*/');
st = [];
en = [];
for i=1:numel(st1)
    en1  = en1(en1>st1(i));
    if ~isempty(en1)
        st  = [st st1(i)];
        en  = [en en1(1)];
        en1 = en1(2:end);
    end
end

str = '';
pen = 1;
for i=1:numel(st)
    str = [str clean_latex(str0(pen:st(i)-1)) str0(st(i)+2:en(i)-1)];
    pen = en(i)+2;
end
str = [str clean_latex(str0(pen:numel(str0)))];
str = remove_cites(str);

%==========================================================================
function str = clean_latex(str)
str  = strrep(str,'$','\$');
%str  = strrep(str,'&','\&');
str  = strrep(str,'^','\^');
%str  = strrep(str,'_','\_');
str  = strrep(str,'#','\#');
%str  = strrep(str,'\','$\\$');
str  = strrep(str,'|','$|$');
%str  = strrep(str,'>','$>$');
%str  = strrep(str,'<','$<$');

function str = remove_cites(str)
% Strip out citations for now
expr = ' *\\cite[^{]*{[^{}]+}';
str  = regexprep(str,expr,'');


%==========================================================================
function bibcstr = get_bib(bibdir)
biblist = dir(fullfile(bibdir,'*.bib'));
bibcstr = {};
for k = 1:numel(biblist)
    n          = spm_file(biblist(k).name,'basename');
    bibcstr{k} = fullfile(bibdir,n);
end
