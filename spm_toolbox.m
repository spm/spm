function spm_toolbox(action,varargin)
% Manage third-party SPM toolboxes
% FORMAT spm_toolbox
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_toolbox.m 4823 2012-08-02 16:23:49Z guillaume $

if ~nargin, action = 'Display'; end
switch lower(action)
    case 'display'
        tbx_display;
    case 'install'
        if isdeployed, return; end
        tbx_install(varargin{:});
    otherwise
        error('Unknown action to perform.');
end

%==========================================================================
% FUNCTION tbx_display
%==========================================================================
function tbx_display

SPMver = spm('Ver');
SPMver = 'SPM8';

%-Initialise window
%--------------------------------------------------------------------------
h = spm_figure('FindWin','SPMtbx');
if isempty(h)
    h = figure(...
        'MenuBar',     'none',...
        'NumberTitle', 'off',...
        'Name',        [SPMver ' Toolboxes'],...
        'Resize',      'off',...
        'Toolbar',     'none',...
        'Tag',         'SPMtbx',...
        'WindowStyle', 'Normal',... 'Modal'
        'Color',       [1 1 1],...
        'Visible',     'off');
    pos = get(h,'Position');
    pos([3 4]) = [350 400];
    set(h,'Position',pos);
else
    pos = get(h,'Position');
end
set(h,'Visible','on');

if ispc, s = '/'; else s = ''; end
localurl = @(f) sprintf(['file://' s strrep(spm('Dir'),'\','/') '/help/' f]);

html_header = sprintf(['<html>\n<head>\n<title>SPM Toolboxes</title>\n',...
    '<link rel="stylesheet" type="text/css" href="%s" />\n',...
    '</head>\n<body>\n'],localurl('spm.css'));
html_footer = sprintf('\n</body>\n</html>');

html = ['<h1 style="text-align:center;">%s Toolboxes</h1>\n' ...
    '<br/><p>Loading description of toolboxes...</p>' ...
    '<p style="margin-top:5em;"><center><img src="%s"></center></p>'];
html = sprintf(html,SPMver,localurl('images/loading.gif'));
html = [html_header html html_footer];
[H, HC] = spm_browser(html,h,[2 2 pos(3)-4 pos(4)-4],'html');
setappdata(h,'browser',HC); set(h,'Resize','on','ResizeFcn',@tbx_resize);

%-Download and parse toolboxes' description
%--------------------------------------------------------------------------
[ext, sts] = urlread('http://www.fil.ion.ucl.ac.uk/spm/ext/ext.xml');
if ~sts
    spm_browser('Download failed.',H);
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
html = sprintf('<h1 style="text-align:center;">%s Toolboxes</h1>\n',SPMver);
ind  = find(ext,'//ext');
nbtbx = 0;
for i=1:numel(ind)
    try, v = attributes(ext,'get',ind(i),'ver'); catch, v = ''; end
    if isempty(v) || ~isempty(strfind(lower(v),lower(SPMver)))
        tbx = convert(ext,ind(i));
        if ismember(tbx.id,tbxinst)
            imgname  = 'images/star2.gif';
            imgtitle = 'Installed';
        else
            imgname  = 'images/star1.gif';
            imgtitle = 'Install';
        end
        html = [html ...
            sprintf(['<div class="tbx">\n<h4><a href="matlab:web(''%s'',''-browser'')">%s</a>' ...
            '<span class="star"><a href="matlab:web(''%s'',''-browser'')"><img src="%s" width="16" height="18" title="%s"/></a></span></h4>\n' ...
            '<p><a href="mailto:%s" title="email %s <%s>">%s</a></p>\n' ...
            '<p>%s</p>\n</div>\n'],...
            tbx.url, tbx.name,...
            tbx.url,localurl(imgname), imgtitle,...
            tbx.email, tbx.author, tbx.email, tbx.author,tbx.summary)];
        nbtbx = nbtbx + 1;
    end
end
html = [html sprintf(['<p>Found %d toolboxes compatible with %s.' ...
    '<span class="star"><a href="matlab:spm_toolbox;"><img src="%s" width="16" height="18" title="Refresh"/></a></span></p>'],...
    nbtbx, SPMver, localurl('images/refresh.gif'))];

html = [html_header html html_footer];

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
H         = getappdata(obj,'browser');
set(H,   'pos',[2 2 figpos(3)-4 figpos(4)-4]);
set(obj, 'Units',old_units);


%==========================================================================
% FUNCTION tbx_install(url,dest)
%==========================================================================
function tbx_install(url,dest)
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
try
    FS = unzip(F,dest);
catch
    error('Could not unzip file "%s".',F);
end
spm_unlink(F);
i = find(~cellfun('isempty', regexpi(FS,'.*README.txt$')));
if ~isempty(i)
    edit(FS{i(1)});
end
