function [H, HC] = spm_browser(url,F,pos,format)
% Display an HTML document within a MATLAB figure
% FORMAT H = spm_browser(url,F,pos,[format])
%
% url    - string containing URL (e.g. 'http://...' or 'file://...')
% F      - figure handle or Tag [Default: Graphics]
% pos    - position within figure in pixel units with the format [left,
%          bottom, width, height]. [Default: full window with 10px margin]
% format - data format {['html'],'wiki'}
%          'wiki' option uses Wikipedia wiki markup:
%          http://en.wikipedia.org/wiki/Wikipedia:Cheatsheet
% 
% H      - handle to the Java component
% HC     - handle to the HG container
%__________________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_browser.m 6409 2015-04-16 16:19:38Z guillaume $

%-Input arguments
%--------------------------------------------------------------------------
if nargin < 1
    url  = 'http://www.fil.ion.ucl.ac.uk/spm/';
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
    format = 'html';
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
    if strcmpi(format,'html') && any(strncmp(url,{'file','http','ftp:'},4))
        H.setCurrentLocation(strrep(url,'\','/'))
    elseif strcmpi(format,'wiki')
        H.setHtmlText(wiki2html(url));
    elseif all(~strncmp(url,{'file','http','ftp:'},4))
        H.setHtmlText(url);
    end
    drawnow
catch
    H  = [];
    HC = [];
end

%==========================================================================
function html = wiki2html(wiki)
% Convert Wiki markup document into HTML

fid = fopen(wiki);
if fid~=-1
    wiki = fscanf(fid,'%c');
    fclose(fid);
else
    html = [sprintf('<html>\n<head>\n<title>wiki2html</title>\n</head>\n'), ...
            sprintf('<body>\n'), ...
            sprintf('<p>Cannot read %s</p>',wiki), ...
            sprintf('\n</body>\n</html>')];
    return;
end

% Section headings
wiki = regexprep(wiki,'======[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?======','<h6>$1</h6>');
wiki = regexprep(wiki,'=====[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?=====','<h5>$1</h5>');
wiki = regexprep(wiki,'====[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?====','<h4>$1</h4>');
wiki = regexprep(wiki,'===[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?===','<h3>$1</h3>');
wiki = regexprep(wiki,'==[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?==','<h2>$1</h2>');
wiki = regexprep(wiki,'=[\s]*?([0-9A-Za-z].[^=\[]*)[\s]*?=','<h1>$1</h1>');

% Horizontal line
wiki = regexprep(wiki,'----','<hr>');

% Inline elements
wiki = regexprep(wiki,'''''''''''([0-9A-Za-z].*)''''''''''',...
    '<strong><em>$1</em></strong>');
wiki = regexprep(wiki,'''''''([0-9A-Za-z].*)''''''','<strong>$1</strong>');
wiki = regexprep(wiki,'''''([0-9A-Za-z].*)''''','<em>$1</em>');

% Paragraphs
wiki = regexprep(wiki,'\n([^#\*=].*)','<p>$1</p>');

% unordered list
wiki = regexprep(wiki,'\*([^*]*)','<li>$1</li>');
wiki = regexprep(wiki,'([^</li>])<li>','$1<ul><li>');
wiki = regexprep(wiki,'<\/li>(?!<li>)','</li></ul>');

% ordered list
wiki = regexprep(wiki,'^#[:]?[#]* (.*)','<li>$1</li>');
wiki = regexprep(wiki,'([^</li>][>]?[\n])<li>','$1<ol><li>');
wiki = regexprep(wiki,'<\/li>\n(?!<li>)','</li></ol>');

% tables
wiki = regexprep(wiki,'^{\|\n\|-','<table><tr>');
wiki = regexprep(wiki,'^\|-','</tr><tr>');
wiki = regexprep(wiki,'^!\s(.*)','<th>$1</th>');
wiki = regexprep(wiki,'^\|\s(.*)','<td>$1</td>');
wiki = regexprep(wiki,'^\|}','</tr></table>');

html = [sprintf('<html>\n<head>\n<title>wiki2html</title>\n</head>\n'), ...
        sprintf('<body>\n'), ...
        wiki, ...
        sprintf('\n</body>\n</html>')];
