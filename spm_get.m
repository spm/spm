function varargout = spm_get(varargin)
% User interface : filename selection
% FORMAT P = spm_get(n,Filter,Prompt,NewWDir,CmdLine);
% n          - 0 - returns with empty P
%            - Positive integer - prompts for n filepaths
%            - +Inf - prompts for filepaths until "Done" pressed
%            - Negative integer - prompts for n directories
%            - -Inf - prompts for directories until "Done" pressed
% Filter     - Filename filter {'*' is prepended} 
%            - numeric 0 for previous filtering
% Prompt     - Prompt string
% NewWDir    - New working directory
% CmdLine    - CmdLine override switch, 0 for GUI, 1 for Command line.
%
% P          - string matrix of full pathnames {one string per row}
%              returned as a cellstr if iscellstr(Prompt)
%_______________________________________________________________________
%
% spm_get allows interactive selection of filepaths/directory names
% from disk.
%
% if global CMDLINE exists and is true, the command line is used. The
% CmdLine parameter overrides this choice.
%
% Otherwise, selection of files is via GUI. The select files window,
% enables the interactive directory display and selection of
% filenames.  Directories are displated in red in the left listing,
% files in the right listing. File lists are filtered according to the
% current filter string, using the usual UNIX conventions (*,?,etc.).
%
% Arrows at the right of the listings pane allow the lists to be
% scrolled up and down. The bar button anove the "up" arrow scrolls to
% the top of the list.
%
% Selecting files
% ---------------
% Clicking on an item (black) adds it to the list of those selected,
% the item is numbered and changes color (dark blue). The last numbered item
% can be unselected by clicking it again. If a specific number of
% items have been requested (n), then only this number can be
% selected. The item window is scrolled up and down with the buttons
% provided. The small button above the ^ button scrolls to the top of the
% listing.
%
% The directory window is editable. Pull down menus of previous
% directories and subdirectories of the current directory are
% maintained. The red "pwd" button changes directory to the current
% MatLab working directory, the "home" button to the users home
% directory. (The bold red list of subdirectories that used to appear
% alongside the files in the item window have had to be withdrawn due
% to a problem with MatLab5.)
%
% All files in the current window can be selected *in the order in which
% they appear* with the "All" button, and the "Reset" button clears the
% current list. The "Edit" button enables interactive editing of the
% filename list, and the "Keybd" button invokes a command line
% interface.
%
% Click "Done" when finished. As a short cut, selecting items
% using the right mouse button is equivalent to selecting the item and
% pressing "Done".
%
% Summary views of directories
% ----------------------------
% Directories with over 48 entries (after filtering) are shown
% summarised, with a single wildcarded entry for blocks of similar
% filenames shown alongside a small white number indicating the number
% of files in the block. Such blocks can be expanded by clicking
% "Extend" (middle mouse button) on the item, or by clicking the number
% of files text, whence the whole directory is shown with the
% appropriate filter. The entire directory can be shown by clicking the
% white "Summary View" text at the top of the item window.
%
% Filter string
% ---------------------
% The filter string for filenames is built as ['*',Filter] if Filter
% doesn't contain a wildcard '*', and appears in a single editable
% window. Subsequent calls to spm_get callbacks with the same filter
% string do not change the value of the filter string, allowing user
% edits of the filer string to persist. Clicking the Filter button
% resets the filter to the original filter string passed by the
% program. Clicking the bar at the right-hand end of the filter window
% sets the filter string to '*', resulting in all files/directories
% being displayed, in summary view if there are sufficient files. (Bars
% above and below the filter window (if available) implement a simple
% "memory" for filter strings.  The top bar stores the current filter
% string, the lower bar recalls it.)
%
% Selecting Directories
% ---------------------
% The same interface is occasionally used for selecting directories, in
% which case the list of selectable items contains a list of
% subdirectories of the current directory. Use the current directory
% window and PullDown menu's to navigate to the required directories
% parent, and select the chosen subdirectory from the list presented.
% The Filter doesn't apply to the directory list.
%
% Keyboard Accelerators
% ---------------------
% Keyboard accelerators are available for item selection, editing the
% current directory and filter windows, and for some of the interface
% buttons:
%
% - Item selection - Items are highlighted (in italic) using ^K & ^J to
% move Up and Down the item list. The highlighted item is (de)selected
% using ^S.
%
% - CurrentDir & Filter editing - Initially, all typed characters, act
% on the CurrentDirectory widget, characters are appended, delete &
% backspace delete the rightmost character. ^E (Edit) toggles between
% the CurrentDir and Filter editable widgets. Return changes to the
% specifiec directory (CurrentDir widget) or applys the current filter
% (Fitler widget). ^U deletes the whole entry, ^W deletes the last
% pathname component, ^I resets the string to it's previous value.
%
% - General Accelerators -
%        Space / ^D   - Done
%        ^P           - pwd
%        ^B           - Keybd
%        Esc          - Reset
%
% - Tab - In addition, the usual tabbing between GUI objects is
% provided by the window interface.
%
% Programmers notes
% -----------------
% Pass positive n for filepaths, negative n for directories.  +Inf/-Inf
% for selection of filepaths/directories until "Done" is pressed.
%
% The 'Select Files Window', created by the first call to spm_get, is
% 'Tag'ged 'SelFileWin'. This window is hidden at the end of the spm_get
% transaction, and is used by subsequent calls. This saves time, and
% also permits the storage of previously visited directories from
% transaction to transaction.
%
% CallBacks are handled as embedded functions within spm_get.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%
% Callback and setup function for spm_get.
% FORMAT varargout = spm_get(varargin)
%           - A multi function function!
%
% FORMAT F = spm_get('CreateFig',LastDirs,Filter)
% Sets up the file selection window as next avaliable figure.
% LastDirs - String matrix for directory history
% Filter   - Filename filter
% F        - Figure used
%
% FORMAT spm_help('FigKeyPressFcn',ch,F)
% Callback for handling keyboard accelerators
% ch       - Character [default get(F,'CurrentCharacter')]
% F        - Figure used [default gcbf]
%
% FORMAT F = spm_get('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.
% F        - Figure used.
% Vis      - Figure visibility (string, 'on'|'off'|'close')
% n        - Number of filepaths/directories to obtain
% Prompt   - Prompt string.
% Filter   - Filename filter.
% WDir     - [optional] New Working directory
% F        - Figure used
%
% FORMAT spm_get('StatusLine',nP,n,F)
% Updates FileSelWin status line
% nP       - #Items already selected
% n        - #Items to select
% F        - Figure used [defaults to 'SelFileWin' 'Tag'ged figure]
%
% FORMAT spm_get('cd',NewDir,bNoDir)
% Callback function for change directory objects
% NewDir   - New Working directory (if not held in current object)
% bNoDir   - Supresses listing of new directory if true
%
% FORMAT spm_get('dir',WDir,Filter,NoComp)
% Lists current working directory and attatches callbacks to contents
% WDir     - Working directory
%            Defaults (unspecified or empty) to UserData of WDir Tagged object
% Filter   - Filename Filter
%            Overrides and sets the Filter object
%            Defaults (unspecified or empty) to UserData of Filter object
% NoComp   - If specified, forces full directory listing - no compression
%
% [sS,I] = FORMAT spm_get('StrSort',S)
% Simple string sort routine for string matrices
% S        - String matrix to be sorted
% sS       - Sorted string matrix
% I        - Sort index: sS = S(I)
%
% FORMAT fS = spm_get('strfliplr',S)
% fliplr for string matrices, ignores trailing padding spaces
% S        - String matrix to be flipped
% fS       - left-right flipped version of S, with last non-space characters
%            of S in first column
%
% FORMAT [FSpecs,FnamePos]=spm_get('FileSummary',Fnames,Cend,Filter,len)
% File summary routine
% Fnames   - String matrix of filenames (all unique)
% Cend     - CommonEND - 'front', 'end' or 'both'
%            Specifys which end to take as common for file summary
%            Defaults to 'both', which gives 'front' precedence over 'back'
% Filter   - Filename filter - suffix filter is appended to 'front' items
% len      - Length of section to define common files
%            Defaults to 3
%
% FORMAT spm_get('Add',h)
% Adds current object (of object h) to selection
% h        - Object handle of an item
%
% FORMAT spm_get('Delete',h)
% Deletes current object from the selection, if it was the last chosen
% h        - Handle of a deletable item (Defaults to gcbo)
%
% FORMAT spm_get('Reset',F)
% Resets the selection and redisplays the directory
% F        - Figure used [defaults to 'SelFileWin' 'Tag'ged figure]
%
% FORMAT spm_get('All',F)
% Adds all remaining unselected items to the selection
% F        - Figure used [defaults to 'SelFileWin' 'Tag'ged figure]
%
% FORMAT P = spm_get('CmdLine',n,Prompt,P,WDir)
% Uses command line to get items
% n        - # items to get
% Prompt   - Prompt string
% P        - Currently selected items
% WDir     - Working directory to use (defaults to pwd)
%
% FORMAT P = spm_get('GUI2CmdLine')
% GUI gateway to command line input.
% P        - Currently selected items
%
% FORMAT spm_get('Edit')
% Invokes window for editing currently selected items
%
% FORMAT spm_get('EditDone',OK)
% Processes the results of the edit.
% OK       - 'OK' or 'CANCEL' - Edit status.
%
% FORMAT spm_get('Done')
% Callback for "Done" button.
%
% FORMAT cpath = spm_get('CPath',path,cwd)
% function to canonicalise paths: Prepends cwd to relative paths, processes
% '..' & '.' directories embedded in path.
% path     - string matrix (or cell array of strings) containing path names
% cwd      - current working directory [defaut '.']
% cpath    - conditioned paths, in same format as input path argument
%
%-----------------------------------------------------------------------
% SUBFUNCTIONS:
%
% FORMAT so = sf_deblank(si)
% Double tailed deblank - strips blanks from both ends of string matrix
%
%_______________________________________________________________________
% Andrew Holmes


%-Parameters
%=======================================================================
%-Default to unlimited file get when no arguments
if nargin == 0, Action=+Inf'; else, Action = varargin{1}; end
PJump = 1;		%-Allow pointer jumping?

if ~ischar(Action)
%=======================================================================
% P = spm_get(n,Filter,Prompt,NewWDir,CmdLine)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<5 CmdLine=[]; else CmdLine=varargin{5}; end
if nargin<4 NewWDir=''; else NewWDir=varargin{4}; end
if nargin<3 Prompt='Select files...'; else Prompt=varargin{3}; end
if nargin<2 | isempty(varargin{2}), Filter=0; else Filter=varargin{2}; end
n = Action;

if (n==0), varargout={[]}; return, end

if isempty(CmdLine)
	global CMDLINE
	if ~isempty(CMDLINE), CmdLine = CMDLINE; else, CmdLine=0; end
end

%-NB: spm_get callbacks use Filter=0 for previous filtering
if ischar(Filter)
	if ~any(Filter=='*'), Filter = ['*',Filter]; end
else
	Filter=0;
end

%-Computation
%=======================================================================

if CmdLine

	%-Command line version, if global CMDLINE is true, & not overridden
	%---------------------------------------------------------------
	P = spm_get('CmdLine',n,char(Prompt),[],NewWDir);

else

	%-GUI version
	%---------------------------------------------------------------
	%-Set up FileSelWin
	[F,cF] = spm_get('Initialise','on',n,char(Prompt),Filter,NewWDir);

	%-Jump cursor to SelFileWin
	if PJump
		PLoc = get(0,'PointerLocation');
		FRec = get(F,'Position');
		set(0,'PointerLocation',[FRec(1)+FRec(3)/2, FRec(2)+FRec(2)/2])
	end
	
	%-Wait until filenames have been selected
	hDone = findobj(F,'Tag','Done');
	waitfor(hDone,'UserData')
	
	%-Exit with error if SelFileWin deleted
	if ~ishandle(hDone), error('SelFileWin deleted!'), end

	%-Get P & exit status
	status = get(hDone,'UserData');
	P = get(findobj(F,'Tag','P'),'UserData');

	%-Reset and hide SelFileWin
	spm_get('Initialise','off');

	%-Return focus to previous figure (if any)
	set(0,'CurrentFigure',cF)

	%-Jump cursor back to previous location
	if PJump, set(0,'PointerLocation',PLoc), end

	%-Exit with error if reset (status=-1)
	if status == -1, error('SPM-GUI reset: spm_get bailing out!'), end
end

%-Convert P to cellstr if Prompt is one
%-----------------------------------------------------------------------
if iscellstr(Prompt), if isempty(P), P={}; else, P = cellstr(P); end, end

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_get : ',char(Prompt),' :'],char(P)); end

varargout = {P};
return
end

%=======================================================================
% - Callbacks & Utility embedded functions
%=======================================================================

switch lower(Action), case 'createfig'
%=======================================================================
% F = spm_get('CreateFig',LastDirs,Filter)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<3 Filter=[]; else Filter=varargin{3}; end
if nargin<2 | isempty(varargin{2}), LastDirs=pwd; else LastDirs=varargin{2}; end
LastDirs=strvcat(LastDirs,getenv('HOME'));
if (exist('spm.m')==2), LastDirs=strvcat(LastDirs,spm('Dir')); end
	
%-Save current figure
cF = get(0,'CurrentFigure');

%-Create window, compute scaling for screen size
%-----------------------------------------------------------------------
S  = get(0,'ScreenSize');
WS = [S(3)/1152 S(4)/900 S(3)/1152 S(4)/900];
F  = figure('IntegerHandle','off',...
	'Tag','SelFileWin',...
	'Name',[spm('GetUser'),' - SPMget'],'NumberTitle','off',...
	'Position',[S(3)/2-400/2,S(4)/2-395/2,400,395].*WS,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'Units','Pixels',...
	'MenuBar','none',...
	'DefaultTextInterpreter','none',...
	'DefaultTextFontName','Helvetica',...
	'DefaultTextFontSize',2*round(12*min(WS)/2),...
	'DefaultUicontrolFontSize',2*round(12*min(spm('GetWinScale'))/2),...
	'DefaultUicontrolInterruptible','on',...
	'Visible','off');
varargout = {F};

%-Reset CurrentFigure to previous figure
set(0,'CurrentFigure',cF)


%-User control objects with callbacks
%-----------------------------------------------------------------------

uicontrol(F,'Style','Frame','Tag','P','UserData',[],...
	'Position',[001 271 400 124].*WS);

% uicontrol(F,'Style','Frame',...
% 	'BackgroundColor',[1 1 1]*.8,...
% 	'Position',[010 355 355 030].*WS);

uicontrol(F,'Style','Text','Tag','Prompt','UserData',[],...
	'String','<Prompt not set yet>',...
	'FontName','Times',...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',2*round(16*min(WS)/2),...
	'ForegroundColor','k',...
	'HorizontalAlignment','Center',...
	'Position',[010 370 380 022].*WS);

uicontrol(F,'Style','PopUp','Tag','LastDirsPopup',...
	'String',strvcat('Previous Directories...',LastDirs),...
	'ToolTipString','cd to previously visited directories',...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'UserData',size(LastDirs,1),...
	'Callback','spm_get(''cd'')',...
	'Position',[010 347 335 022].*WS);

uicontrol(F,'Style','PushButton','Tag','CDpwd',...
	'String','pwd',...
	'ToolTipString','change to current Matlab working directory',...
	'ForegroundColor','r',...
	'UserData',0,...
	'CallBack','spm_get(''cd'',pwd)',...
	'Position',[345 347 045 022].*WS);

WDir=deblank(LastDirs(1,:));
uicontrol(F,'Style','Edit','String',WDir,...
	'ToolTipString','current directory - edit (& press RETURN) to cd',...
	'Tag','WDir','UserData',WDir,...
	'ForegroundColor','r','BackgroundColor',[.8,.8,1],...
	'HorizontalAlignment','Left',...
	'CallBack','spm_get(''cd'')',...
	'Position',[010 325 380 022].*WS);

uicontrol(F,'Style','PopUp','Tag','SubDirsPopup',...
	'ToolTipString','cd to subdirectories of this directory',...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'String','SubDirectories...',...
	'Callback','spm_get(''cd'')',...
	'Position',[010 303 335 022].*WS);

uicontrol(F,'Style','PushButton','Tag','CDhome',...
	'String','home',...
	'ToolTipString','cd to your home directory',...
	'ForegroundColor','r',...
	'CallBack',['spm_get(''cd'',''',deblank(LastDirs(2,:)),''')'],...
	'Position',[345 303 045 022].*WS);

% uicontrol(F,'Style','Pushbutton','Tag','FilterStore',...
% 	'String','',...
% 	'CallBack',[...
% 		'set(findobj(gcbf,''Tag'',''FilterButton''),''UserData'',',...
% 		    'get(findobj(gcbf,''Tag'',''Filter''),''String'')),'],...
% 	'Position',[010 300 140 005].*WS);

uicontrol(F,'Style','Pushbutton','Tag','FilterButton',...
	'String','Filter:',...
	'ToolTipString','reset filter to ''*''',...
	'CallBack','spm_get(''dir'',[],''*'')',...
	'UserData','*',...
	'Position',[010 278 040 022].*WS);

% uicontrol(F,'Style','Pushbutton','Tag','FilterRecall',...
% 	'String','',...
% 	'CallBack',[...
% 		'spm_get(''dir'',[],get(findobj(gcbf,''Tag'',',...
% 		'''FilterButton''),''UserData''))'],...
% 	'Position',[010 273 140 005].*WS);

uicontrol(F,'Style','Pushbutton','Tag','Filter*',...
	'String','',...
	'ToolTipString','reset filter to initial value',...
	'CallBack',[...
		'spm_get(''dir'',[],get(findobj(gcbf,''Tag'',',...
			'''Filter''),''UserData''))'],...
	'Position',[145 278 005 022].*WS);

uicontrol(F,'Style','Edit','Tag','Filter',...
	'String','*',...
	'ToolTipString','UNIX style file filter string - edit to change',...
	'UserData','*',...
	'BackgroundColor',[.8,.8,1],...
	'Callback','spm_get(''dir'')',...
	'Position',[050 278 095 022].*WS);

uicontrol(F,'Style','Pushbutton','String','All',...
	'ToolTipString','select all items in this directory',...
	'Callback','spm_get(''All'',gcbf)',...
	'Position',[155 278 030 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Edit',...
	'ToolTipString','edit list of items',...
	'ForegroundColor','b',...
	'Callback','spm_get(''Edit'')',...
	'Position',[190 278 040 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Keybd',...
	'ToolTipString','switch to keyboard input of items',...
	'ForegroundColor','k',...
	'Callback','spm_get(''GUI2CmdLine'')',...
	'Position',[235 278 050 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Reset',...
	'ToolTipString','reset selection',...
	'ForegroundColor','r',...
	'Callback','spm_get(''Reset'',gcbf)',...
	'Position',[290 278 050 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Done',...
	'ToolTipString','done - press when completed selection of items',...
	'ForegroundColor','m',...
	'Tag','Done','UserData',1,...
	'Callback','spm_get(''Done'')',...
	'Interruptible','off','BusyAction','Cancel',...
	'Position',[345 278 045 022].*WS);

uicontrol(F,'Style','Pushbutton','String','',...
	'ToolTipString','scroll to top',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''UserData''))'],...
	'Position',[370 265 020 005].*WS);

uicontrol(F,'Style','Pushbutton','String','/\',...
	'ToolTipString','scroll up',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''Position'')-[0 0.5 0 0])'],...
	'Position',[370 245 020 020].*WS);

uicontrol(F,'Style','Pushbutton','String','\/',...
	'ToolTipString','scroll down',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''Position'')+[0 0.5 0 0])'],...
	'Position',[370 040 020 020].*WS);

uicontrol(F,'Style','Frame','Tag','StatusArea',...
	'Position',[001 000 400 030].*WS);

if exist('spm_help.m')==2
	uicontrol(F,'Style','Pushbutton','String','?',...
		'ToolTipString','help on using spm_get',...
		'ForegroundColor','g',...
		'Callback','spm_help(''spm_get.m'')',...
		'Position',[370 005 020 020].*WS);
end

uicontrol(F,'Style','Text','Tag','StatusLine',...
	'String','<Not set yet>',...
	'FontAngle','Italic',...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[010 005 355 020].*WS);

set(F,'KeyPressFcn',...
	'spm_get(''FigKeyPressFcn'',get(gcbf,''CurrentCharacter''),gcbf)')


case 'figkeypressfcn'
%=======================================================================
% spm_help('FigKeyPressFcn',ch,F)
if nargin<3, F = gcbf; else, F = varargin{3}; end
if nargin<2, ch=get(F,'CurrentCharacter'); else, ch=varargin{2}; end

%-Empty ch - shift/control/&c.
if isempty(ch), return, end

%-Keyboard accelerators for item selection
%-----------------------------------------------------------------------
if any(abs(ch)==[10,11,19]) % ^J,^K,^S
	H = flipud(findobj(get(get(F,'CurrentAxes'),'Children'),...
		'Flat','HandleVisibility','on'));
	if isempty(H), return, end
	h = max([0,findobj(H,'FontWeight','bold')]);
	hPos = max([0;find(H==h)]);
	if abs(ch)==19 & hPos
		%- ^S - Select
		switch get(h,'Tag')
			case 'IName', spm_get('Add',h)
			case 'SelIName', spm_get('Delete',h)
			case 'DirName', spm_get('cd',get(h,'UserData'))
		end
		return
	end
	nhPos = max(1,min(hPos-2*(abs(ch)-10.5),length(H)));
	if hPos, set(h,'FontWeight','normal'), end
	set(H(nhPos),'FontWeight','bold')
	return
end

%-General Accelerators
%-----------------------------------------------------------------------
if any(abs(ch)==[32,4,16,2,27]) % Space,^D,^P,^B,Esc
	switch abs(ch)
	case {32,4}	%- Space/^D - Done
		spm_get('Done')
	case 16		%- ^P - Change to MatLab pwd
		spm_get('cd',pwd)
	case 2		%- ^B - Keyboard
		spm_get('GUI2CmdLine')
	case 27		%- Esc - Reset
		spm_get('Reset')
	end
end

%-Accelerators for WDir & Filter Edit widgets
%-----------------------------------------------------------------------
%-Edit which widget? WDir or Filter?
hCDpwd    = findobj(F,'Tag','CDpwd');
EditWDir  = get(hCDpwd,'UserData');
if EditWDir, h = findobj(F,'Tag','WDir');
	else, h = findobj(F,'Tag','Filter'); end
str = get(h,'String');

if abs(ch)==9		%- ^I - reset string
	str = get(h,'UserData');
elseif abs(ch)==5	%- ^E - switch between Edit widgets, reset strings
	set(hCDpwd,'UserData',~EditWDir)
	str = get(h,'UserData');
elseif abs(ch)==23	%- ^W - word delete
	tmp = max(find(str=='/'));
	if ~isempty(tmp), str(tmp:end)=[]; end
elseif abs(ch)==13	%- Return - Goto topic / apply filter
	if EditWDir, spm_get('cd',str)
		else, spm_get('dir'), end
	return
elseif any(abs(ch)==[32:126])
	str = [str, ch];
elseif abs(ch)==21	%- ^U - kill
	str = '';
elseif any(abs(ch)==[8,127])	%-BackSpace or Delete
	if length(str), str(length(str))=''; end
else	%-Other character
	return
end
set(h,'String',str)


case 'initialise'
%=======================================================================
% [F,cF] = spm_get('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.

if nargin<6 WDir=''; else WDir=varargin{6}; end
if nargin<5 Filter=0; else Filter=varargin{5}; end
if nargin<4 Prompt='Select files...'; else Prompt=varargin{4}; end
if nargin<3 n=Inf; else n=varargin{3}; end
if nargin<2 Vis='on'; else Vis=varargin{2}; end


%-Recover SelFileWin figure number
F  = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

cF = get(0,'CurrentFigure');				%-Save current figure

switch lower(Vis), case 'close'
	close(F)
	varargout = {[],cF};
	return
case 'reset'
	varargout = {F,cF};				%-Return figure handles
	if isempty(F), return, end
	set(F,'Visible','off')				%-Make window Invisible
	set(findobj(F,'Tag','Done'),'UserData',-1)	%-Set Done UserData to -1
	delete(get(F,'CurrentAxes'))			%-delete 'dir' axes
	drawnow
	return
case 'off'
	if isempty(F), varargout={[],cF}; return, end
	varargout = {F,cF};				%-Return figure handles
	set(F,'Visible','off')				%-Make window Invisible
	set(findobj(F,'Tag','Done'),'UserData',1)	%-Set Done UserData to 1
	delete(get(F,'CurrentAxes'))			%-delete 'dir' axes
	drawnow
	return
case 'on'
	if isempty(F), 	F=spm_get('CreateFig'); end	%-Create SelFileWin
	varargout = {F,cF};				%-Return figure handles
	set(findobj(F,'Tag','Done'),'UserData',0)	%-Init. Done UserData
	delete(get(F,'CurrentAxes'))			%-delete 'dir' axes
	set(findobj(F,'Tag','P'),'UserData','')		%-clear P
	set(findobj(F,'Tag','Prompt'),'String',Prompt,'UserData',n)
							%-Reset prompt, & n
	spm_get('StatusLine',0,n,F)			%-Reset StatusLine
	set(findobj(F,'Tag','CDpwd'),'UserData',1)	%-KeyPressFcn to edit
							% WDir widget initially
	%-Filter string: Only set if different to previously passed one
	% (held as UserData of Filter Edit uicontrol)
	if ischar(Filter)
	  if ~strcmp(Filter,get(findobj(F,'Tag','Filter'),'UserData'))
	    set(findobj(F,'Tag','Filter'),'String',Filter,'UserData',Filter)
	  end
	end

	%-Change to new working directory, if specified, else display current
	if ~isempty(WDir)
		spm_get('cd',WDir)			%-Change directory
	else
		spm_get('dir')				%-Show current dir.
	end

	%-Popup figure, retaining CurrentFigure
	figure(F)					%-PopUp figure
	set(0,'CurrentFigure',cF)			%-Return to prev. figure

otherwise
	error('Unrecognised ''Vis'' option in spm_get(''Initialise'',...')
end



case 'statusline'
%=======================================================================
% spm_get('StatusLine',nP,n,F)
if nargin<4, F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
	else, F = varargin{4}; end
if nargin<3, n=get(findobj(F,'Tag','Prompt'),'UserData');else,n=varargin{3}; end
if nargin<2, nP=size(get(findobj(F,'Tag','P'),'UserData'),1);else,nP=varargin{2}; end

str=['Selected ',int2str(nP)];
if finite(n), str=[str,'/',int2str(abs(n))]; end
if (n>=0)
	str = [str,' file'];
	if nP~=1, str=[str,'s']; end
else
	str = [str,' director'];
	if nP==1, str=[str,'y']; else, str=[str,'ies']; end
end

if isinf(n)
	str = [str,', press "Done" when finished.'];
elseif nP==abs(n)
	str = [str,', press "Done".'];
else
	str = [str,'.'];
end

set(findobj(F,'Tag','StatusLine'),'String',str)


case 'cd'
%=======================================================================
% spm_get('cd',NewDir,bNoDirList)

%-Condition arguments
%-----------------------------------------------------------------------
if (nargin<3), bDirList=1; else, bDirList = ~varargin{3}; end
if (nargin<2), NewDir=''; else, NewDir=varargin{2}; end

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(F,'Pointer','Watch')
WDir=get(findobj(F,'Tag','WDir'),'UserData');

if nargin<2
	%-Current object contains NewDir, either UIEdit, UIPopMenu or text
	NewDir=get(gcbo,'String');
	if strcmp(get(gcbo,'Style'),'popupmenu')
		val = get(gcbo,'Value');
		if val==1 set(F,'Pointer','Arrow'), return, end
		NewDir=NewDir(val,:);
	end
else
	NewDir=varargin{2};
end


%-Condition directory path
%-----------------------------------------------------------------------
NewDir = spm_get('CPath',NewDir,WDir);


%-Changing directory, delete current directory listing (Done again in
% Action=='dir', but neater to do here in advance on this occasion)
%-----------------------------------------------------------------------
delete(get(F,'CurrentAxes'))


%-Set up LastDirs
%-----------------------------------------------------------------------
%-Recover LastDirs & size of fixed part from popup menu object
LastDirs  = get(findobj(F,'Tag','LastDirsPopup'),'String'); LastDirs(1,:)=[];
Fs        = get(findobj(F,'Tag','LastDirsPopup'),'UserData');

%-Add new directory
LastDirs  = strvcat(NewDir,LastDirs);

%-Extract fixed and variable parts of LastDirs
FLastDirs = LastDirs(end-Fs+1:end,:);
LastDirs(end-Fs+1:end,:)=[];

%-Delete any replications of NewDir within the variable part of LastDirs
if size(LastDirs,1)>1
	IDRows = 1+find(all( LastDirs(2:end,:) == ...
		repmat(LastDirs(1,:),size(LastDirs,1)-1,1) ,2 ) );
	LastDirs(IDRows,:)=[];
end

%-Delete NewDir from top of LastDirs if present in fixed part of LastDirs
if any(all(FLastDirs==repmat(LastDirs(1,:),Fs,1),2)), LastDirs(1,:)=[]; end

%-Lose directories from the bottom of variable LastDirs if oversize
Vs = 10; %-Maximum allowable number of variable LastDirs entries
if size(LastDirs,1)>Vs, LastDirs(Vs+1:end,:)=[]; end

%-Put LastDirs and FLastDirs back together
LastDirs = [LastDirs; FLastDirs];


%-Write LastDirs to LastDirsPopup, write new directory to WDir Edit widget
%-----------------------------------------------------------------------
set(findobj(F,'Tag','LastDirsPopup'),'Value',1,...
	'String',strvcat('Previous Directories...',LastDirs))
set(findobj(F,'Tag','WDir'),'String',NewDir,'UserData',NewDir)


%-Display new directory, unless supressed
%-----------------------------------------------------------------------
if bDirList
	spm_get('dir')
else
	set(findobj(F,'Tag','SubDirsPopup'),'Value',1,...
		'String','SubDirectories...')
	set(F,'Pointer','Arrow')
	drawnow
end


case 'dir'
%=======================================================================
% spm_get('dir',WDir,Filter,NoComp)
% Creates list of text objects with an associated
% ButtonDownFcn functions that call spm_get callbacks for processing
%
%-WDir - Directory, defaults to UserData of 'WDir' tagged object in
% SelFileWin figure

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(F,'Pointer','Watch')
delete(get(F,'CurrentAxes')), drawnow

%-Condition parameters and setup variables
%-----------------------------------------------------------------------
if nargin<4, NoComp=0; else, NoComp=1; end
if nargin<3, Filter=''; else, Filter=varargin{3}; end
if nargin<2 | isempty(varargin{2})
	WDir = get(findobj(F,'Tag','WDir'),'UserData');
else
	WDir = deblank(varargin{2});
end

n      = get(findobj(F,'Tag','Prompt'),'UserData');

%-Sort out Filter
%-----------------------------------------------------------------------
if isempty(Filter)
	Filter = get(findobj(F,'Tag','Filter'),'String');
	if isempty(Filter), Filter='*'; end
end
set(findobj(F,'Tag','Filter'),'String',Filter)

%-Write current directory to WDir widget, if new directory specified
%-----------------------------------------------------------------------
if nargin>1
	set(findobj(F,'Tag','WDir'),'String',WDir,'UserData',WDir)
end

%-Set up axis
%-----------------------------------------------------------------------
Rec = [0.02 0.086 0.90 0.58];
hAxes = axes('Parent',F,...
	'Position',Rec,...
	'DefaultTextInterpreter','none',...
	'DefaultTextHandleVisibility','off',...
	'Units','Points',...
	'Visible','off',...
	'UserData',Rec,...
	'Tag','hAxes');
y     = floor(get(hAxes,'Position'));
y0    = y(3);
dy    = 22;
set(hAxes,'Ylim',[0,y0])


%-List current directory
%-----------------------------------------------------------------------
[Files,Dirs] = spm_list_files(WDir,Filter);
if isempty(Dirs)
	text(0,y0,'Permission denied, or non-existent directory',...
		'Parent',hAxes,'FontWeight','bold','Color','r');
	set(F,'Pointer','Arrow')
	return
else	%-Exclude "." directory from Dirs, exclude ".." if in /
	%-Exclude ".." from DirItems (used when selecting directories)
	Dirs(Dirs(:,1)=='.',:)=[];
	DirItems = strvcat('.',Dirs);
	if ~strcmp(WDir,'/'), Dirs=strvcat('..',Dirs); end
end

%-Create list of directories in pulldown menu
%-----------------------------------------------------------------------
set(findobj(F,'Tag','SubDirsPopup'),'Value',1,...
	'String',strvcat('SubDirectories...',Dirs))

%-Create list of directories with appropriate 'ButtonDownFcn': -
% NB: Gave MatLab Bus errors (ML5.0.0.4064, Solaris2.5.1)
%-----------------------------------------------------------------------
y     = y0-dy;
for i = 1:size(Dirs,1)
	text(0,y,deblank(Dirs(i,:)),...
		'Parent',hAxes,...
		'Tag','DirName',...
		'FontWeight','bold','Color','r',...
		'UserData',deblank(Dirs(i,:)),...
		'ButtonDownFcn','spm_get(''cd'',get(gcbo,''UserData''))',...
		'HandleVisibility','off');
	y = y - dy;
end

%-Files or directories (n<0)
%-----------------------------------------------------------------------
if n>0	%-Select file(s) - omit ".*" dot files
	if ~isempty(Files), Files(Files(:,1)=='.',:)=[]; end
	Items = Files;
else	%-Select directory/ies
	Items = DirItems;
end

%-Compressed summary view, or full view?
%-----------------------------------------------------------------------
if ~NoComp & size(Items,1)>32
	%-Use a compressed summary view
	if Filter(end)=='*'
		[IName,ItemPos] = spm_get('FileSummary',Items);
	else
		[IName,ItemPos] = spm_get('FileSummary',Items,'front',Filter);
	end
	text(0.5,y0,'Summary View - Click to expand',...
		'Parent',hAxes,...
		'Color','w',...
		'FontSize',10,...
		'HorizontalAlignment','Center',...
		'FontAngle','Italic',...
		'ButtonDownFcn','spm_get(''dir'',[],[],1)',...
		'HandleVisibility','off');
else
	%-Use a standard view, each item with it's own representation
	IName = Items; ItemPos = [1:size(Items,1)]';
end

%-Create list of Items with appropriate 'ButtonDownFcn': -
%-----------------------------------------------------------------------
y     = y0-dy;
for i = 1:size(IName,1)
	cIName   = IName(i,IName(i,:)~=' ');
	cItemPos = ItemPos(i,ItemPos(i,:)>0);
	nItems   = length(cItemPos);
	cItems   = Items(cItemPos,:);
	%-Next line strips off redundant space at end of string matrix
	% (Don't have to worry about spaces within strings 'cos filenames)
	cItems(:,all(cItems==' ',1))=[];
	text(0.35,y,cIName,...
		'Parent',hAxes,...
		'Tag','IName',...
		'UserData',cItems,...
		'Color','k',...
		'ButtonDownFcn','spm_get(''Add'')',...
		'HandleVisibility','on')
	if nItems>1;
		text(0.34,y-3,int2str(sum(ItemPos(i,:)>0)),...
			'Parent',hAxes,...
			'Color','w',...
			'FontSize',9,...
			'FontAngle','Italic',...
			'HorizontalAlignment','right',...
			'UserData',cIName,...
			'ButtonDownFcn',...
			'spm_get(''dir'',[],get(gcbo,''UserData''),1)',...
			'HandleVisibility','off')
	end
	y = y - dy;
end
set(F,'Pointer','Arrow')


case 'strsort'
%=======================================================================
% [sS,I]=spm_get('StrSort',S)
% Utility string sorting routine for string matrices
if nargin<2, error('Sort what?'), else, S = varargin{2}; end

if size(S,2)<=3
	a        = min(abs(S(:)'));
	b        = max(abs(S(:)'))+1;
	i        = (b-a).^[size(S,2)-1:-1:0];
	[null,I] = sort(abs(S)*i');
	sS       = S(I,:);
else
	lS            = size(S,2);
	hS            = S(:,1:3);
	tS            = S;
	tS(:,1:3)     = [];
	[null,tI]     = spm_get('StrSort',tS);
	[null,hI]     = spm_get('StrSort',hS(tI,:));
	I             = tI(hI);
	sS            = S(I,:);	
end
varargout = {sS,I};


case 'strfliplr'
%=======================================================================
% fS = spm_get('strfliplr',S)
if nargin<2, error('Flip what?'), else, S = varargin{2}; end

%-FlipLR the Fnames matrix, but watch out for trailing blanks/spaces
[nS,Sl] = size(S);
fS = ones(nS,Sl);
for i = 1:nS
	c = max(find(S(i,:)~=' '));
	fS(i,:) = S(i,[c:-1:1,Sl:-1:c+1]);
end
varargout = {fS};


case 'filesummary'
%=======================================================================
% [FSpecs,FnamePos]=spm_get('FileSummary',Fnames,Cend,Filter,len)
if nargin<5, len = 3; else len = varargin{5}; end
if nargin<4 | isempty(varargin{4}), Filter = '*'; else Filter = varargin{4}; end
if nargin<3, Cend='both'; else, Cend = varargin{3}; end
if nargin<2 | isempty(varargin{2}), error('Specify Fnames to summarise')
else Fnames = varargin{2}; end
%-Implicit assumption in this code is that no two files have the same name

if strcmp(Cend,'front')
	if size(Fnames,1) == 1	%-Only one filename!
		varargout = {deblank(Fnames),1};
		return
	end
	WildCard    = Filter(max(find(Filter=='*')):end);
	nFnames     = size(Fnames,1);
	[sFnames,I] = spm_get('StrSort',Fnames);
	tmp         = diff(abs(sFnames));
	i           = [0,find(any(tmp(:,1:len)')),nFnames];
	nFgroups    = length(i)-1;
	FSpecs      = '';
	FnamePos    = [];
	for Fgroup  = 1:nFgroups
		cI       = [i(Fgroup)+1:i(Fgroup+1)];
		Cfnames  = sFnames(cI,:);
		nCfnames = size(Cfnames,1);
		if nCfnames == 1
			tmp    = max(find(Cfnames~=' '));
			FSpecs = strvcat(FSpecs,Cfnames(1,1:tmp));
		elseif nCfnames == 2
			tmp    = min(find([diff(abs(Cfnames)),1]))-1;
			FSpecs = strvcat(FSpecs,[Cfnames(1,1:tmp),WildCard]);
		else
			tmp    = min(find([any(diff(abs(Cfnames))),1]))-1;
			FSpecs = strvcat(FSpecs,[Cfnames(1,1:tmp),WildCard]);
		end
		tI       = I(cI);
		tmp	 = max([nCfnames,size(FnamePos,2)]);
		FnamePos = [FnamePos,...
				zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
				tI', zeros(1,tmp-length(tI))];
	end
elseif strcmp(Cend,'end')
	fFnames           = spm_get('strfliplr',Fnames);
	[FSpecs,FnamePos] = spm_get('FileSummary',fFnames,'front');
	FSpecs            = spm_get('strfliplr',FSpecs);
elseif strcmp(Cend,'both')
	[hSpecs,hI] = spm_get('FileSummary',Fnames,'front');
	FSpecs      = '';
	FnamePos    = [];
	for i = 1:size(hSpecs,1)
	    chI          = hI(i,hI(i,:)>0);
	    if length(chI) == 1
		%-Single file
		FSpecs   = strvcat(FSpecs,hSpecs(i,hSpecs(i,:)~=' '));
		tmp      = max([1,size(FnamePos,2)]);
		FnamePos = [FnamePos,...
			zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
			chI, zeros(1,tmp-1)];
	    else
		%-Multiple files - process by common endings
		Cfnames     = Fnames(chI,:);
		%-Chop off common beginning
		Cfnames(:,1:find(hSpecs(i,:)=='*')-1) = [];
		%-Watch out for blank lines (e.g for files abc & abcde)
		bI = find(all(Cfnames==' ',2));
		if ~isempty(bI)
		    % There'll be one blank line at most
		    FSpecs   = strvcat(FSpecs,...
			hSpecs(i,1:find(hSpecs(i,:)=='*')-1) );
		    tmp      = max([1,size(FnamePos,2)]);
		    FnamePos = [FnamePos,...
			zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
			chI(bI), zeros(1,tmp-1)];
		    Cfnames(bI,:) = [];
		    chI(bI)       = [];
		end
		%-Even if we've deleted a blank line, there'll be some left
		[tSpecs,tI] = spm_get('FileSummary',Cfnames,'end');
		for j = 1:size(tSpecs,1)
		    %-Sort in next line keeps ordering by headers
		    cI       = chI(sort(tI(j,tI(j,:)>0)));
		    if length(cI) == 1
			%-Unique, don't bother with '*'s
			cFSpec = [hSpecs(i,1:find(hSpecs(i,:)=='*')-1),...
				tSpecs(j,1:max(find(tSpecs(j,:)~=' ')) ) ];
		    else
		        %-Don't double '*'s
		        cFSpec = [hSpecs(i,hSpecs(i,:)~=' '),...
				tSpecs(j,2:max(find(tSpecs(j,:)~=' ')) ) ];
		    end
		    FSpecs   = strvcat(FSpecs,cFSpec);
		    tmp      = max([length(cI),size(FnamePos,2)]);
		    FnamePos = [FnamePos,...
			zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
			cI, zeros(1,tmp-length(cI))];
		end
	    end
	end
else
	error('Invalid Cend specifier')
end
varargout = {FSpecs,FnamePos};


case 'add'
%=======================================================================
% spm_get('Add',H)
% Add filename to current list
%-H - [Optional] (vector of) handles to file/dir name objects
%

%-Recover variables from holding objects
%-----------------------------------------------------------------------
F       = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
SelType = get(F,'SelectionType');
n       = get(findobj(F,'Tag','Prompt'),'UserData');
P       = get(findobj(F,'Tag','P'),'UserData');
WDir    = get(findobj(F,'Tag','WDir'),'UserData');

if (nargin==2)
	%-Multiple handles passed by 'All' CallBack
	H=varargin{2};
else
	%-'Add' called by item callback - get handle & check selection
	H=gcbo;
	if strcmp(SelType,'extend')
		IName = get(H,'String');
		if any(IName=='*')
			spm_get('dir','',IName,1)
			return
		else
			return
		end
	elseif strcmp(SelType,'alt')
		spm_get('Add',H)
		spm_get('Done')
		return
	end
end

%-Add items to P (until P would be overfull)
%-----------------------------------------------------------------------
for h = H'
	IName  = get(h,'String');
	Items  = get(h,'UserData');
	nItems = size(Items,1);
	nP     = size(P,1);

	%-Don't add any more if #items has been specified,
	% and P would be overfull
	if finite(n), if (nP+nItems>abs(n)), break, end, end

	%-Add items to P
	FPath  = [repmat([WDir,'/'],nItems,1),Items];
	P      = strvcat(P,FPath);
	if nItems==1
		tmp = [int2str(nP+1),' :',IName];
	else
		tmp = [int2str(nP+1),'-',int2str(nP+nItems),' :',IName];
	end
	set(h,'String',tmp,'Color','b',...
		'Tag','SelIName',...
		'ButtonDownFcn','spm_get(''Delete'')')
end

%-Chop off any redundant trailing spaces
%-----------------------------------------------------------------------
% (Don't have to worry about spaces within strings 'cos filenames)
if ~isempty(P), P(:,all(P==' ',1))=[]; else, P=''; end

%-Return new P to holding object
%-----------------------------------------------------------------------
set(findobj(F,'Tag','P'),'UserData',P);

%-Update StatusLine
%-----------------------------------------------------------------------
spm_get('StatusLine',size(P,1),n,F)


case 'delete'
%=======================================================================
% spm_get('Delete',h)
F      = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
if nargin<2, h=gcbo; else, h=varargin{2}; end
WDir   = get(findobj(F,'Tag','WDir'),'UserData');
P      = get(findobj(F,'Tag','P'),'UserData');
n      = get(findobj(F,'Tag','Prompt'),'UserData');
IName  = get(h,'String');
IName(1:find(IName==':')) = [];
Items  = get(h,'UserData');
nItems = size(Items,1);
FPath  = [repmat([WDir,'/'],nItems,1),Items];
nP     = size(P,1);

%-Compare Items with end of P (Can only delete from the end)
tmp    = strvcat(P(nP-(nItems-1):nP,:),FPath);
if all(all(tmp(1:nItems,:)==tmp(nItems+1:2*nItems,:)))
	P(nP-(nItems-1):nP,:)=[];
	if ~isempty(P), P(:,all(P==' ',1))=[]; else, P=''; end
	set(findobj(F,'Tag','P'),'UserData',P);
	set(h,'String',IName,'Color','k',...
		'Tag','IName',...
		'ButtonDownFcn','spm_get(''Add'')')
	spm_get('StatusLine',size(P,1),n,F)
end


case 'reset'
%=======================================================================
% spm_get('Reset',F)
if nargin<2, F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
	else, F = varargin{2}; end
set(findobj(F,'Tag','P'),'UserData',[]);
spm_get('dir')
spm_get('StatusLine')


case 'all'
%=======================================================================
% spm_get('All',F)
if nargin<2, F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
	else, F = varargin{2}; end
set(F,'Pointer','Watch')
H = flipud(findobj(get(get(F,'CurrentAxes'),'Children'),'Flat','Tag','IName'));
spm_get('Add',H);
set(F,'Pointer','Arrow')


case 'cmdline'
%=======================================================================
% P = spm_get('CmdLine',n,Prompt,P,WDir,AllowEnd)
if nargin<6, AllowEnd=0; else AllowEnd=varargin{6}; end
if nargin<5, WDir=''; else WDir=varargin{5}; end
if isempty(WDir), WDir=pwd; else, WDir=deblank(WDir); end
if nargin<4, P=[]; else P=varargin{4}; end
if nargin<3, Prompt='Select files...'; else Prompt=varargin{3}; end
if nargin<2, n=Inf; else n=varargin{2}; end

if n==0, return, end

nP=size(P,1);

fprintf('\n'), fprintf('%c',repmat('=',1,70)), fprintf('\n')
fprintf('\t%s\n',Prompt)
fprintf('%c',repmat('=',1,70)), fprintf('\n')
fprintf('Current working directory: %s\n',WDir)
fprintf('(Prepended to relative paths)\n')
if nP>0
	fprintf('\nAlready selected:\n')
	for pP = [1:size(P,1)], fprintf('\t%s\n',P(pP,:)), end
end
fprintf('\nEnter paths :')
AllowEnd = AllowEnd | isinf(n);
Tstr = 'END';
if AllowEnd, fprintf(' (Type "END" to terminate input.)\n')
	else fprintf(' to %d items:\n',abs(n)-nP), end

Done=0;
while ~Done
	%-Get input
	str = []; while isempty(str)
		str=sf_deblank(input(sprintf('  %3d  : ',nP+1),'s'));end
	%-Prepend WDir to incomplete pathnames
	if (str(1)~='/')&(~strcmp(str,Tstr)), str = [WDir,'/',str]; end
	if (n>0),OK=exist(str)==2;else,OK=~unix(['test -d ',str]);end
	while ~( OK | (strcmp(str,Tstr) & AllowEnd) )
		if strcmp(str,Tstr)
			fprintf('%c\tSelect %d files!',7,abs(n))
			else, fprintf('%c\t%s doesn''t exist!',7,str), end
		str=[]; while isempty(str)
			str=sf_deblank(input(sprintf('  %3d  : ',nP+1),'s'));end
		if (str(1)~='/')&(~strcmp(str,Tstr)), str=[WDir,'/',str]; end
		if (n>0),OK=exist(str)==2;else,OK=~unix(['test -d ',str]);end
	end
	if ~strcmp(str,Tstr), P=strvcat(P,str); end
	nP = nP+1;
	Done = (strcmp(str,Tstr) | (nP==abs(n)));
end

varargout = {P};


case 'gui2cmdline'
%=======================================================================
% P = spm_get('GUI2CmdLine')

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
Prompt = get(findobj(F,'Tag','Prompt'),'String');
n      = get(findobj(F,'Tag','Prompt'),'UserData');
WDir   = get(findobj(F,'Tag','WDir'),'UserData');
P      = get(findobj(F,'Tag','P'),'UserData');

WS = get(F,'Position'); WS = [WS(3)/400,WS(4)/395,WS(3)/400,WS(4)/395];
Handles = [];
h = uicontrol(F,'Style','Frame',...
		'Position',[000 000 400 395].*WS);
Handles = [Handles, h];

h = uicontrol(F,'Style','Text','String',Prompt,...
	'FontName','Times',...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',2*round(16*min(WS)/2),...
	'ForegroundColor','k',...
	'HorizontalAlignment','Center',...
	'Position',[010 370 380 022].*WS);
Handles = [Handles, h];

h = uicontrol(F,'Style','Frame',...
	'Position',[001 001 400 030].*WS);
Handles = [Handles, h];

h = uicontrol(F,'Style','Text',...
	'String','Use command window...',...
	'FontAngle','Italic',...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[020 005 360 020].*WS);
Handles = [Handles, h];

PNew = spm_get('CmdLine',n,Prompt,P,WDir,1);

if ~strcmp(P(:),PNew(:))
	set(findobj(F,'Tag','P'),'UserData',PNew);
	spm_get('dir')
	spm_get('StatusLine',size(PNew,1),n,F)
end

delete(Handles);


case 'edit'
%=======================================================================
% spm_get('Edit')

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
n = get(findobj(F,'Tag','Prompt'),'UserData');

set(F,'Units','Pixels');
WS = get(F,'Position'); WS = [WS(3)/400,WS(4)/395,WS(3)/400,WS(4)/395];

EditHandles = [];

h = uicontrol(F,'Style','Frame',...
	'Tag','EditHandles',...
	'Position',[001 001 400 395].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Text',...
	'String',get(findobj(F,'Tag','Prompt'),'String'),...
	'FontName','Times',...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',2*round(16*min(WS)/2),...
	'ForegroundColor','k',...
	'HorizontalAlignment','Center',...
	'Position',[010 370 380 022].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Text',...
	'String','Edit the filename window...',...
	'ForegroundColor','w',...
	'HorizontalAlignment','Left',...
	'Position',[010 347 210 022].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Pushbutton','String','Cancel',...
	'ForegroundColor','r',...
	'Callback','spm_get(''EditDone'',''Cancel'')',...
	'Position',[220 347 070 022].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Pushbutton','String','OK',...
	'ForegroundColor','m',...
	'Callback','spm_get(''EditDone'',''OK'')',...
	'Position',[310 347 070 022].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Edit',...
	'String',get(findobj(F,'Tag','P'),'UserData'),...
	'Tag','EditWindow','UserData',0,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','b','BackgroundColor',[.8,.8,1],...
	'Max',2,...
	'Callback','set(gco,''UserData'',1)',...
	'Position',[010 040 380 300].*WS);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Frame',...
	'Position',[001 001 400 030].*WS);
EditHandles = [EditHandles, h];

str=[]; if finite(n), str=['Use only top ',int2str(abs(n)),' lines.  ']; end
h = uicontrol(F,'Style','Text',...
	'String',[str,'(Blank lines will be removed)'],...
	'FontAngle','Italic',...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[020 005 360 020].*WS);
EditHandles = [EditHandles, h];

set(findobj(F,'Tag','EditHandles'),'UserData',EditHandles)


case 'editdone'
%=======================================================================
% spm_get('EditDone',OK)
if nargin<2, OK=0; else, OK=strcmp(varargin{2},'OK'); end
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

EditHandles = get(findobj(F,'Tag','EditHandles'),'UserData');
set(EditHandles,'Visible','off')
h_EditWindow = findobj(F,'Tag','EditWindow');

if OK & get(h_EditWindow,'UserData')
	%-Only update if OK was pressed, and there was an edit.
	P = get(h_EditWindow,'String');
	n = get(findobj(F,'Tag','Prompt'),'UserData');
	%-Condition P, by removing blank lines and truncating to n items.
	if ~isempty(P), P = P(~all(P'==0),:); end
	if (finite(n))&(~isempty(P)), P=P(1:min(size(P,1),abs(n)),:); end
	%-Write P to 'P' 'Tag'ged object
	set(findobj(F,'Tag','P'),'UserData',P);
	%-Redo directory listing
	spm_get('dir')
	%-Redo status line
	spm_get('StatusLine',size(P,1),n,F)
end
delete(EditHandles), drawnow


case 'done'
%=======================================================================
% spm_get('Done')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
n = get(findobj(F,'Tag','Prompt'),'UserData');

%-If #files was specified, and not enough have been selected, then return
if finite(n)
	P = get(findobj(F,'Tag','P'),'UserData');
	if (size(P,1)<abs(n)), fprintf('%c',7), return, end
end

%-Done, set Done UserData tag for handling
set(findobj(F,'Tag','Done'),'UserData',1)


case 'cpath'
%=======================================================================
% cpath = spm_get('cpath',path,cwd)
if nargin<3, cwd=pwd; else, cwd=spm_get('CPath',varargin{3}); end
if nargin<2, error('insufficient arguments'), end
cpath = cellstr(varargin{2});

for i=1:prod(size(cpath))

	%-Prepend cwd to relative pathnames
	%---------------------------------------------------------------
	if isempty(cpath{i}) | cpath{i}(1)~='/'
		if strcmp(cwd,'/'), cpath{i}=['/',cpath{i}];
			else cpath{i}=[cwd,'/',cpath{i}]; end
	end
	
	%-Sort out stationary relative pathnames './' & '/.'
	%---------------------------------------------------------------
	%-Remove midpath '/./'
	cpath{i}=strrep(cpath{i},'/./','/');
	while length(cpath{i})>=2 & length(findstr(cpath{i},'//'))
		cpath{i}=strrep(cpath{i},'//','/'); end
	
	%-Remove '.' from trailing '/.'
	if length(cpath{i})>=2 & strcmp(cpath{i}(end-1:end),'/.')
		cpath{i}(end)=''; end
	
	%-Remove trailing '/'
	if length(cpath{i})>=2 & cpath{i}(end)=='/'
		cpath{i}(end)=''; end
	
	
	%-Sort out relative pathnames '/..'
	%---------------------------------------------------------------
	%-Canonicalise paths starting '/..' to '/'
	if length(cpath{i})>=3 & strcmp(cpath{i}(1:3),'/..')
		cpath{i}(2:3)=''; end
	
	%-Process midpath '/../'
	t0 = 0;
	while length(cpath{i})>=4 & max([0,findstr(cpath{i},'/../')])>t0
		t2 = t0+min(findstr(cpath{i}(t0+1:end),'/../'));
		t1 = t0+max([0,find(cpath{i}(t0+1:t2-1)=='/')]);
		if t1==t0 | strcmp(cpath{i}(t1:t2),'/../')
			t0=t2;
		else
			cpath{i}(t1:t2+2)='';
		end
	end
	
	%-Process trailing '/..'
	if length(cpath{i})>3 & strcmp(cpath{i}(end-2:end),'/..')
		t2 = length(cpath{i})-2;
		t1 = t0+max([0,find(cpath{i}(t0+1:t2-1)=='/')]);
		if t1~=t0 & ~strcmp(cpath{i}(t1:t2),'/../')
			cpath{i}(t1:t2+2)='';
		end
	end
	
	%-Final canonicalisations
	%---------------------------------------------------------------
	% %-Remove leading './'
	% if length(cpath{i})>2 & strcmp(cpath{i}(1:2),'./')
	% 	cpath{i}(1:2)=''; end
	
	%-Canonicalise leading './../' to '../'
	if length(cpath{i})>5 & strcmp(cpath{i}(1:5),'./../')
		cpath{i}(1:2)=''; end
	%-Canonicalise './..' to '..'
	if strcmp(cpath{i},'./..'), cpath{i}='..'; end

end

if ischar(varargin{2}), varargout={char(cpath)}; else, varargout={cpath}; end


otherwise
%=======================================================================
error('Illegal Action string')


%=======================================================================
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function so = sf_deblank(si)
%=======================================================================
if nargin==0, error('insufficient arguments'), end
if isempty(si), so=''; return, end
if ~ischar(si), error('input must be a string'), end

b = any((si~=' ' & si~=0),1);
if all(~b), so=''; return, end
so = si(:,max(1,min(find(b))):min(max(find(b)),length(b)))
