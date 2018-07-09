function varargout = spm_toolbox(action,varargin)
% Manage third-party SPM toolboxes
% FORMAT spm_toolbox
% FORMAT spm_toolbox install <tbxname>
%
% See https://www.fil.ion.ucl.ac.uk/spm/ext/
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_toolbox.m 7372 2018-07-09 16:50:44Z guillaume $


if ~nargin, action = 'Display'; end
switch lower(action)
    case 'display'
        tbx_display;
    case 'install'
        if isdeployed, return; end
        tbx_install(varargin{:});
    case {'update', 'uninstall', 'list', 'describe', 'build'}
        error('Not implemented yet.');
    case 'form'
        tbx_form(varargin{:});
    case 'url'
        varargout = {'http://www.fil.ion.ucl.ac.uk/spm/ext/ext.xml'};
    case 'spmver'
        varargout = {spm('Ver')};
    otherwise
        error('Unknown action to perform.');
end

%==========================================================================
% FUNCTION tbx_display
%==========================================================================
function tbx_display

SPMver = spm_toolbox('SPMver');

if ispc, s = '/'; else s = ''; end
localurl = @(f) sprintf(['file://' s strrep(spm('Dir'),'\','/') '/help/' f]);

%-Initialise HTML template
%--------------------------------------------------------------------------
tpl = spm_file_template(fullfile(spm('Dir'),'help'),'keep');
tpl = tpl.file('TPL_TBX','spm_toolbox.tpl');
tpl = tpl.block('TPL_TBX','loading','load');
tpl = tpl.block('TPL_TBX','listing','list');
tpl = tpl.block('listing','toolbox','tbx');
tpl = tpl.var('SPM',SPMver);
tpl = tpl.var('SPM_CSS',localurl('spm.css'));

%-Initialise window
%--------------------------------------------------------------------------
h   = tbx_figure([SPMver ' Toolboxes']);
H   = getappdata(h,'Hbrowser');

tpl = tpl.var('IMG_LOADING',localurl('images/loading.gif'));
tpl = tpl.parse('load','loading',0);
tpl = tpl.var('list','');
tpl = tpl.parse('OUT','TPL_TBX');
html = get(tpl,'OUT');

spm_browser(html,H);

%-Download and parse toolboxes' description
%--------------------------------------------------------------------------
[ext, sts] = urlread(spm_toolbox('url'));
if ~sts
    spm_browser('Server not found. Cannot list SPM toolboxes.',H);
    return
end
spm_browser(strrep(html,'Loading','Parsing'),H);
ext = xmltree(ext);

%-Detect already installed toolboxes
%--------------------------------------------------------------------------
tbxdir  = spm_get_defaults('tbx.dir');
tbxinst = {};
for i=1:numel(tbxdir)
    d   = dir(tbxdir{i});
    d   = {d([d.isdir]).name};
    d(strncmp('.',d,1)) = [];
    tbxinst = [tbxinst d];
end

%-Select relevant toolboxes
%--------------------------------------------------------------------------
ind  = find(ext,'//ext');
nbtbx = 0;
tpl = tpl.var('load','');
tpl = tpl.var('tbx','');
for i=1:numel(ind)
    try, v = attributes(ext,'get',ind(i),'ver'); catch, v = ''; end
    if isempty(v) || ~isempty(strfind(lower(v),lower(SPMver)))
        tbx = convert(ext,ind(i));
        if ismember(tbx.id,tbxinst)
            imgname  = 'images/star2.gif';
            imgtitle = 'Installed';
            % matlab:web(''%s'',''-browser'')
        else
            imgname  = 'images/star1.gif';
            imgtitle = 'Install';
            % matlab:spm_toolbox(''install'',''%s'')
        end
        nbtbx = nbtbx + 1;
        tpl = tpl.var('TBX_URL',tbx.url);
        tpl = tpl.var('TBX_NAME',tbx.name);
        tpl = tpl.var('TBX_ID',tbx.id);
        tpl = tpl.var('IMG_STAR',localurl(imgname));
        tpl = tpl.var('IMG_STAR_TITLE',imgtitle);
        tpl = tpl.var('TBX_EMAIL',tbx.email);
        tpl = tpl.var('TBX_AUTHOR',tbx.author);
        tpl = tpl.var('TBX_SUMMARY',tbx.summary);
        tpl = tpl.parse('tbx','toolbox',1);
    end
end
tpl = tpl.var('TBX_NB',int2str(nbtbx));
tpl = tpl.var('IMG_REFRESH',localurl('images/refresh.gif'));
tpl = tpl.parse('list','listing',0);
tpl = tpl.parse('OUT','TPL_TBX');
html = get(tpl,'OUT');

%-Display
%--------------------------------------------------------------------------
spm_browser(html,H);


%==========================================================================
% FUNCTION tbx_resize(obj,evt,varargin)
%==========================================================================
function tbx_resize(obj,evt,varargin)
old_units = get(obj,'Units');
set(obj, 'Units','pixels');
figpos    = get(obj,'Position');
H         = getappdata(obj,'HCbrowser');
set(H,   'pos',[2 2 figpos(3)-4 figpos(4)-4]);
set(obj, 'Units',old_units);


%==========================================================================
% FUNCTION tbx_install(tbx,dest)
%==========================================================================
function tbx_install(tbx,dest)

if ~nargin, spm_toolbox; return; end

SPMver = spm_toolbox('SPMver');
SPMurl = spm_toolbox('url');

%-Get toolbox details
%--------------------------------------------------------------------------
[ext, sts] = urlread(SPMurl);
if ~sts, error('Cannot access "%s".',SPMurl); end
ext  = xmltree(ext);
iExt = find(ext,'/extensions/ext');
vId  = get(ext,children(ext,find(ext,'/extensions/ext/id')),'value');
iTbx = find(ismember(lower(vId),lower(tbx)));
if isempty(iTbx), error('Cannot find toolbox "%s".',tbx); end
iExt = iExt(iTbx(1));
try, v = attributes(ext,'get',iExt,'ver'); catch, v = ''; end
if ~isempty(v) && isempty(strfind(lower(v),lower(SPMver)))
    warning('Toolbox "%s" is not knowingly compatible with %s.',tbx,SPMver);
end
ext  = convert(ext,iExt);
try
    url = ext.download.url;
    if isempty(url), error('Empty URL.'); end
catch
    str = ext.url;
    if spm_platform('desktop'), str = sprintf('<a href="%s">%s</a>',str,str); end
    fprintf(['Toolbox "%s" is only available through website:\n' ...
        '  ' str '\n'], tbx);
    if ~spm('CmdLine'), web(ext.url,'-browser'); end
    %spm_toolbox('form',tbx);
    return;
end

%-Check destination directory
%--------------------------------------------------------------------------
if nargin == 1
    tbxdir = spm_get_defaults('tbx.dir');
    dest   = tbxdir{end};
end
[sts, attrb] = fileattrib(dest);
if ~sts, error('"%s"\n%s',dest,attrb); end
if ~attrb.UserWrite
    error('No write access to "%s".\nMaybe use "%s" instead.',...
        dest, strrep(userpath,pathsep,''));
end

%-Download toolbox archive
%--------------------------------------------------------------------------
tmpfile = [tempname(dest) '.zip'];
try
    F = urlwrite(url,tmpfile);
catch
    l = lasterror;
    switch l.identifier
        case 'MATLAB:urlwrite:ConnectionFailed'
            error('Could not access URL "%s".',url);
        case 'MATLAB:urlwrite:InvalidOutputLocation'
            error('Could not create output file "%s".',tmpfile);
        otherwise
            rethrow(l);
    end
end

%-Unzip archive in destination folder
%--------------------------------------------------------------------------
try
    FS = unzip(F,dest);
catch
    spm_unlink(F);
    error('Error when unpackig toolbox archive');
end

%-Delete toolbox archive
%--------------------------------------------------------------------------
spm_unlink(F);

%-Display README.txt or Contents.m if present
%--------------------------------------------------------------------------
i = find(~cellfun('isempty', regexpi(FS,'.*README.txt$')));
if isempty(i)
    i = find(~cellfun('isempty', regexpi(FS,'.*Contents.m$')));
end
if ~isempty(i)
    if ~spm('CmdLine')
        open(FS{i(1)});
    else
        type(FS{i(1)});
    end
end

fprintf('Toolbox "%s" installed in "%s".\n',tbx,dest);

%==========================================================================
% FUNCTION tbx_form(tbxname)
%==========================================================================
function tbx_form(tbxname)

h = tbx_figure(sprintf('Toolbox %s',tbxname));
H = getappdata(h,'Hbrowser');

html = sprintf([...
'<p style="font-weight:bold;">The authors of "%s" request you to fill in this form.</p>\n' ...
'<form method="post" action="matlab:spm_toolbox;">\n' ...
'    <div style="">\n' ...
'        <div>Name</div>\n' ...
'        <div><input type="text" name="name" value="" size="50" /></div>\n' ...
'        \n' ...
'        <div>Email</div>\n' ...
'        <div><input type="text" name="email" value="" size="50" /></div>\n' ...
'        \n' ...
'        <div>Institution</div>\n' ...
'        <div><input type="text" name="institution" value="" size="50" /></div>\n' ...
'        \n' ...
'        <div>Country</div>\n' ...
'        <div><input type="text" name="country" value="" size="50" /></div>\n' ...
'        \n' ...
'        <div>Comments</div>\n' ...
'        <div><textarea name="comments" rows="5" cols="40" value=""></textarea></div>\n' ...
'        \n' ...
'        <div><input type="checkbox" name="licence" value="agree" checked="checked" />\n' ...
'        I accept the terms and conditions of the licence agreement.</div>\n' ...
'        \n' ...
'        <div><input type="hidden" name="SPM" value="%s" /></div>\n' ...
'        \n' ...
'        <div><input value="Proceed to download" type="submit" /></div>\n' ...
'    </div>\n' ...
'</form>\n'],tbxname,spm('Ver'));

spm_browser(html,H);

%==========================================================================
% FUNCTION tbx_figure
%==========================================================================
function h = tbx_figure(name)

if ~nargin, name = 'SPM Toolboxes'; end
h = spm_figure('FindWin','SPMtbx');
if ~isempty(h), set(h,'Name',name); return; end

h = figure(...
    'MenuBar',     'none',...
    'NumberTitle', 'off',...
    'Name',        name,...
    'Resize',      'off',...
    'Toolbar',     'none',...
    'Tag',         'SPMtbx',...
    'WindowStyle', 'Normal',... %'Modal'
    'Color',       [1 1 1],...
    'Visible',     'off');
pos = get(h,'Position');
pos([3 4]) = [350 400];
set(h,'Position',pos);

[H, HC] = spm_browser('<html></html>',h,[2 2 pos(3)-4 pos(4)-4],'html');
setappdata(h,'Hbrowser',H); setappdata(h,'HCbrowser',HC);
set(h,'Resize','on','ResizeFcn',@tbx_resize);

set(h,'Visible','on');
