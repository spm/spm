function spm_browser(url,format)
% Display an HTML document in a web browser
% FORMAT spm_browser(url,[format])
%
% url    - string containing URL (e.g. 'http://...' or 'file://...')
% format - data format {['html'],'md'}
%          'md' option uses Markdown:
%            https://www.wikipedia.org/wiki/Markdown
%            http://showdownjs.com/
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging


%-Input arguments
%--------------------------------------------------------------------------
if nargin < 1
    url  = ['file://' fullfile(spm('Dir'),'help','index.html')];
end

if nargin < 2 || isempty(format)
    if numel(url) > 1 && strcmp(url(end-2:end),'.md')
        format = 'md';
    else
        format = 'html';
    end
end

%-Display
%--------------------------------------------------------------------------
if strcmpi(format,'html') && any(strncmp(url,{'file','http'},4))
    web(strrep(url,'\','/'),'-browser');
else
    if strcmpi(format,'md')
        url = md2html(url);
    end
    if ~isdeployed
        url = ['text://' url];
        web(url,'-notoolbar','-noaddressbox');
    else
        tmp = [tempname '.html'];
        spm_save(tmp,url);
        web(tmp,'-browser');
    end
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
