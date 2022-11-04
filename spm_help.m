function varargout = spm_help(varargin)
% SPM help and manual facilities
% FORMAT spm_help
%
% The "Help" facilities are about software and implementation. The
% underlying mathematics, concepts and operational equations have been (or
% will be) published in the peer reviewed literature and the interested
% user is referred to these sources. An intermediate theoretical exposition
% is given in the SPM course notes. This and other resources are available
% via the SPM Web site.
% Visit https://www.fil.ion.ucl.ac.uk/spm/, or press the "SPMweb" button.
%
%--------------------------------------------------------------------------
%
% `spm_help('Topic')` or `spm_help Topic` displays the help for a
% particular topic.
%
%--------------------------------------------------------------------------
% The SPM package provides help at three levels, the first two being
% available via the SPM graphical help system:
%
% (i)   Manual pages on specific topics.
%       These give an overview of specific components or topics its
%       relation to other components, the inputs and outputs and
%       references to further information.
%
%       Many of the buttons in the help menu window lead to such "man"
%       pages.  These are contained in Markdown files named spm_*.md.
%       These can be viewed on the MATLAB command line with the `help`
%       command, e.g. `help spm_help` prints out this manual file in
%       the MATLAB command window.
%
% (ii)  Help information for each routine within SPM (E.g. This is the).
%       help information for spm_help.m - the help function.)
%       This help information is the help header of the actual MATLAB
%       function, and can be displayed on the command line with the
%       `help` command, e.g. `help spm_help`.
%
% (iii) SPM is (mainly) implemented as MATLAB functions and scripts.
%       These are ASCII files named spm_*.m, which can be viewed in the
%       MATLAB command window with the `type` command, e.g. `type
%       spm_help`, or read in a text editor.
%
%  ---  MATLAB syntax is very similar to standard matrix notation that
%       would be found in much of the literature on matrices. In this
%       sense the SPM routines can be used (with MATLAB) for data
%       analysis, or they can be regarded as the ultimate pseudocode
%       specification of the underlying ideas.
%
% In addition, the MATLAB help system provides keyword searching through
% the H1 lines (the first comment line) of the help entries of *all*
% M-files found on MATLABPATH. This can be used to identify routines from
% keywords. Type `help lookfor` in the MATLAB command window for further
% details.
%__________________________________________________________________________

% Andrew Holmes, Karl Friston
% Copyright (C) 1994-2022 Wellcome Centre for Human Neuroimaging


%__________________________________________________________________________
%
% FORMAT spm_help
% Defaults to spm_help('!Topic').
%
% FORMAT spm_help('Topic')
% Defaults to spm_help('!Topic','Topic').
%
% FORMAT spm_help('!Topic',Topic)
% Topic     - Help topic: Either name of file from which to display help,
%             or name of an internal help topic
%             Defaults to README.md
% Loads file Topic and displays it in the Help window.
%
% FORMAT spm_help('!Disp',Fname,S)
% Fname     - Name of file from which to display help
% S         - [Optional] String vector containing a previously read in
%             contents of file Fname
% Displays the help for the given file in the Help window (creating
% one if required).
%__________________________________________________________________________


%-Condition arguments
%--------------------------------------------------------------------------
%-All actions begin '!' - Other (string) actions are topics
if nargin==0 || isempty(varargin{1})
    spm_help('!Topic');
    return
elseif varargin{1}(1)~='!'
    spm_help('!Topic',varargin{:}); return
end

%-
%--------------------------------------------------------------------------
action = varargin{1};
switch lower(action)

case {'!quit','!drawmenuwin','!shorttopics','!create','!createhelpwin',...
'!createbar','!pulldowntopic','!figkeypressfcn','!clear','!contexthelp'}
%==========================================================================
    warning('Action ''%s'' is deprecated.',action);

case '!disp'
% FORMAT spm_help('!Disp',Fname,S)
%==========================================================================
    if nargin<3, S='';     else S=varargin{3}; end
    if nargin<2, Fname = fullfile(spm('Dir'),'README.md');
    else         Fname = varargin{2}; end
    
    if isempty(S), S = get_content(Fname); end
    spm_browser(S);

case '!topic'
% FORMAT spm_help('!Topic',Topic)
%==========================================================================
    if nargin > 1
        topic = varargin{2};
    else
        topic = fullfile(spm('Dir'),'README.md');
    end
    spm_help('!Disp',topic);
    
otherwise
%==========================================================================
    error('Unknown action string');

end


%==========================================================================
% function url = get_content(topic)
%==========================================================================
function url = get_content(topic)
% Get content to be displayed concerning topic

switch lower(spm_file(topic,'ext'))
    case 'm'
        H = help(topic);
    case {'html','htm'}
        if any(strncmp(topic,{'file','http','ftp:'},4))
            url = topic;
        else
            url = ['file:///' topic];
        end
    case {'txt','man'}
        fid = fopen(topic,'rt');
        if fid == -1, H = '';
        else H = fread(fid,'*char')'; fclose(fid); end
        H = strrep(H,char([10 37]),char(10));
        if numel(H) && H(1)=='%', H(1)=''; end
    case 'md'
        url = topic;
    otherwise
        if exist([topic '.m'],'file')
            url = get_content([topic '.m']);
        else
            H = spm_jobman('help',topic,80);
            [H{2,:}] = deal(sprintf('\n'));
            H = [H{:}];
        end
end
if ~exist('url','var')
    url = ['<h1>' spm_file(topic,'filename') '</h1>',...
        '<pre style="font-size:15;">' escapechar(H) '</pre>'];
end


%==========================================================================
% function f = fullurl(varargin)
%==========================================================================
function f = fullurl(varargin)
% Build URL from parts

f = strrep(fullfile(varargin{:}),'\','/');


%==========================================================================
% function str = escapechar(str)
%==========================================================================
function str = escapechar(str)
% Escape HTML special characters
% See http://www.w3.org/TR/html4/charset.html#h-5.3.2

str = strrep(str,'&','&amp;');
str = strrep(str,'<','&lt;');
str = strrep(str,'>','&gt;');
str = strrep(str,'"','&quot;');


%==========================================================================
% function varargout = spmfcn
%==========================================================================
function varargout = spmfcn
% List filenames of all SPM functions

persistent b f
if ~isempty(b) && ~isempty(f), varargout = {b,f}; return; end

f = cellstr(spm_select('FPListRec',spm('Dir'),'^.*\m$'));
b = cellfun(@(x) spm_file(x,'basename'),f,'UniformOutput',false);
i = cellfun(@(x) ~isempty(strmatch('spm_',x)),b);
f = f(i);
b = b(i);

varargout = {b,f};
