function [H, HC] = spm_browser(url,F,pos,format)
% Display an HTML document within a MATLAB figure
% FORMAT H = spm_browser(url,F,pos,[format])
%
% url    - string containing URL (e.g. 'http://...' or 'file://...')
% F      - figure handle or Tag [Default: Graphics]
% pos    - position within figure in pixel units with the format [left,
%          bottom, width, height]. [Default: full window with 10px margin]
% format - data format {['html'],'md'}
%          'md' option uses Markdown:
%            https://www.wikipedia.org/wiki/Markdown
%            http://showdownjs.com/
% 
% H      - handle to the Java component
% HC     - handle to the HG container
%__________________________________________________________________________
% Copyright (C) 2011-2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_browser.m 7478 2018-11-08 14:51:54Z guillaume $

%-Input arguments
%--------------------------------------------------------------------------
if nargin < 1
    url  = ['file://' fullfile(spm('Dir'),'help','index.html')];
end

if nargin < 2 || isempty(F)
    F    = spm_figure('GetWin','Graphics');
end
try, isUpdate = ismethod(F,'setCurrentLocation'); catch, isUpdate=false; end
if ~isUpdate, F = spm_figure('FindWin',F); end

if (nargin < 3 || isempty(pos)) && ~isUpdate
    u    = get(F,'Units');
    set(F,'Units','pixels');
    pos  = get(F,'Position');
    set(F,'Units',u);
    pos  = [10, 10, pos(3)-20, pos(4)-20];
end

if nargin < 4 || isempty(format)
    if numel(url) > 1 && strcmp(url(end-2:end),'.md')
        format = 'md';
    else
        format = 'html';
    end
end

%-Display
%--------------------------------------------------------------------------
try
    % if usejava('awt') && spm_check_version('matlab','7.4') >= 0
    %-Create HTML browser panel
    %----------------------------------------------------------------------
    if ~isUpdate
        browser = com.mathworks.mlwidgets.html.HTMLBrowserPanel;
        [H, HC] = javacomponent(browser,pos,F);
    else
        H       = F;
        HC      = [];
    end
    
    %-Set content
    %----------------------------------------------------------------------
    if strcmpi(format,'html') && any(strncmp(url,{'file','http'},4))
        H.setCurrentLocation(strrep(url,'\','/'))
    elseif strcmpi(format,'md')
        H.setHtmlText(md2html(url));
    else
        H.setHtmlText(url);
    end
    drawnow;
catch
    H  = [];
    HC = [];
end


%==========================================================================
function html = md2html(md)
% Convert Markdown document into HTML (using Showdown.js)

if exist(md,'file')
	md = fileread(md);
elseif any(strncmp(md,{'file','http'},4))
    md = urlread(md);
end
md = strrep(md,'\','\\');
md = strrep(md,'"','\"');
md = strrep(md,char(13),'');
md = strrep(md,char(10),'\n');

showdownjs = ['file://' fullfile(spm('Dir'),'help','js','showdown.min.js')];

html = {...
    '<!DOCTYPE html>', ....
    '<html>', ....
      '<head>', ....
        '<meta charset="utf-8"/>', ....
        '<title>SPM</title>', ....
        ['<script src="', showdownjs,'"></script>'], ....
      '</head>', ....
      '<body>', ....
        '<div id="display">Loading...</div>', ....
        '<script>', ....
          'var converter = new showdown.Converter();', ....
          'converter.setFlavor("github");', ...
          'converter.setOption("simpleLineBreaks", false);', ...
          ['var text   = "', md,'";'], ....
          'var html    = converter.makeHtml(text);', ....
          'var display = document.getElementById("display");', ....
          'display.innerHTML = html;', ....
        '</script>', ....
      '</body>', ....
    '</html>'};
html = sprintf('%s\n', html{:});
