function varargout = spm_get(varargin)
% User interface : filename selection
% FORMAT P = spm_get(n,Filter,Prompt,NewWDir,CmdLine);
% n          - 0 - returns with empty P
%            - [mn,mx]		: prompts for between mn & mx items
%            - scalar n => [n,n], i.e. must select n items
%            - all(n) positive  : prompts for files
%            - any(n) negative  : prompts for directories
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
% if global CMDLINE exists and is >0, the command line is used. The
% CmdLine parameter overrides this choice. Note that since spm_input
% uses the command line if CMDLINE is true (unless overridden), setting
% CMDLINE to -1 results in the command line being used for questions,
% but the GUI for file selection.
%
% Otherwise, selection of files is via GUI. The select files window,
% enables the interactive directory display and selection of
% filenames.  Directories are displated in red in the left listing,
% files in the right listing. File lists are filtered according to the
% current filter string, using the usual UNIX conventions (*,?,etc.).
%
% Arrows at the right of the listings pane allow the lists to be
% scrolled up and down. The bar button above the "up" arrow scrolls to
% the top of the list.
%
% Selecting files
% ---------------
% Clicking on an item (black) adds it to the list of those selected,
% the item is numbered and changes color (dark blue). The item last
% selected can be unselected by clicking it again. If a specific number
% of items have been requested (n), then this number of files must be
% selected. A status message at the foot of the window explains the
% current status.
%
% Clicking on a directory name (red) changes to that directory. The
% current directory name window can edited to change directory. Also,
% pull down menus of previous directories and subdirectories of the
% current directory are maintained. (On Windows platforms, an
% additional pull down menu of drives is also provided.) The red "pwd"
% button changes directory to the current MatLab working directory, the
% "home" button to the users home directory.
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
% resets the filter to '*', resulting in all files/directories being
% displayed, in summary view if there are sufficient files. Clicking
% the outline of the filter editable widget resets the filter string to
% the original filter string passed by the calling program.
%
% Filter can also be a cell array of the form {'*.img', 'noexpand'}.
% When it is a cell array, the first element is taken as the filter
% to use, and subsequent elements are passed as the third argument to
% spm_list_files_expand.m.
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
% buttons: These work providing no GUI widget has the focus (click in
% the window background to remove focus from any GUI widgets).
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
% The 'Select Files Window', created by the first call to spm_get, is
% 'Tag'ged 'SelFileWin'. This window is hidden at the end of the spm_get
% transaction, and is used by subsequent calls. This saves time, and
% also permits the storage of previously visited directories from
% transaction to transaction.
%
% CallBacks are handled as embedded functions within spm_get, format
% specifications for these embedded functions are given in the main
% body of the code.
%
%_______________________________________________________________________
% @(#)spm_get.m	2.45 Andrew Holmes (X-platform stuff with Matthew Brett) 04/05/17

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
% n        - Items struct: .fi (get files), .mn (min #items), .mx (max #items)
% Prompt   - Prompt string.
% Filter   - Filename filter.
% WDir     - [optional] New Working directory
% F        - Figure used
%
% FORMAT spm_get('StatusLine',nP,n,F)
% Updates FileSelWin status line
% nP       - #Items already selected
% n        - Items struct: .fi (get files), .mn (min #items), .mx (max #items)
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
%            Specifies which end to take as common for file summary
%            Defaults to 'both', which gives 'front' precedence over 'back'
% Filter   - Filename filter - suffix filter is appended to 'front' items
% len      - Length of section to define common files
%            Defaults to 3
% FSpecs   - String matrix of file patterns, with n rows
% FnamePos - Matrix of n rows.
%            Each row contains the indicies of all Fnames matching the pattern
%            given in the corresponding row of FSpecs. Padded with zeros.
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
% n        - Items struct: .fi (get files), .mn (min #items), .mx (max #items)
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
% FORMAT str = spm_get('DrivesPullDownStr')
% returns string for Drives pull down menu on Windows platforms
%
% FORMAT [P,dir] = spm_get('Files',dir,fil)
% Extended gateway to spm_list_files: Returns full pathnames & canonicalises dir
% dir      - directory within which list files (a string)
%            parameter is canonicalised by spm_get('CPath',dir)
%            [defaults to '.', the current directory]
% fil      - file filter string E.g. 'sn*.img' [default '*']
% P        - string matrix of files in directory
%            if nargout<2 (i.e. if dir is not returned), then these are
%            full pathnames
% dir      - full pathname of directory listed (after canonicalisation)
%
%-----------------------------------------------------------------------
% SUBFUNCTIONS:
%
% FORMAT so = ddeblank(si)
% Double tailed deblank - strips blanks from both ends of string matrix
%
% FORMAT rootf = isroot(path)
% Tests filepath in string path: Returns true if path is root directory
%
% FORMAT absf = isabspath(path)
% Tests filepath in string path: Returns true if path is absolute
%_______________________________________________________________________
% Andrew Holmes (X-platform stuff with Matthew Brett)

% Options to pass to spm_list_files_expand
%=======================================================================
persistent ListOpts
if isempty(ListOpts), ListOpts = {''}; end;

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
CmdLine = spm('CmdLine',CmdLine);
if nargin<4 NewWDir=''; else NewWDir=varargin{4}; end
if nargin<3 Prompt='Select files...'; else Prompt=varargin{3}; end
if nargin<2 | isempty(varargin{2}),
	Filter=0;
else
	Filter   = varargin{2};
	ListOpts = {''};
	if any(Action(:)<0),
		ListOpts = {'noexpand'};
	end;
	if iscell(Filter),
		ListOpts = Filter(2:end);
		Filter   = Filter{1};
	end;
end

%-Parse number of items parameter (Passed as "Action" parameter)
%-----------------------------------------------------------------------
n  = struct(	'fi',	~any(Action(:)<0),...		%-Get files flag
		'mn',	min(abs(Action(:))),...		%-Min # items
		'mx',	max(abs(Action(:)))	);	%-Max # items
if isinf(n.mn), n.mn=0; end

if n.mx==0
	if iscellstr(Prompt), varargout={{}}; else, varargout={''}; end
	return
end

%-NB: spm_get callbacks use Filter=0 for previous filtering
%-----------------------------------------------------------------------
if ischar(Filter)
	if ~any(Filter=='*'), Filter = ['*',Filter]; end
else
	Filter=0;
end

%-Computation
%=======================================================================

if CmdLine>0

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
	if ~ishandle(hDone), error('SPM file selector window was quit!'), end

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
	spm_log([mfilename,' : ',char(Prompt),' :'],P); end

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
if nargin<3 Filter=[];
else
	Filter=varargin{3};
end
if nargin<2 | isempty(varargin{2}), LastDirs=pwd; else LastDirs=varargin{2}; end
LastDirs=strvcat(LastDirs,getenv('HOME'));
if (exist('spm.m')==2), LastDirs=strvcat(LastDirs,spm('Dir')); end

%-Save current figure
cF = get(0,'CurrentFigure');

%-Create window, compute scaling for screen size
%-----------------------------------------------------------------------
WS = spm('WinScale');				%-Window scaling factors
FS = spm('FontSizes');				%-Scaled font sizes
PF = spm_platform('fonts');			%-Font names (for this platform)
S0 = get(0,'ScreenSize');			%-Screen size

F  = figure('IntegerHandle','off',...
	'Tag','SelFileWin',...
	'Name',sprintf('%s%s: SPMget',spm('ver'),spm('GetUser',' (%s)')),...
	'NumberTitle','off',...
	'Position',[S0(3),S0(4),0,0]/2 + [-400/2,-395/2,400,395].*WS,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'Units','Pixels',...
	'MenuBar','none',...
	'DefaultTextInterpreter','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
 	'DefaultAxesFontName',PF.helvetica,...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(10),...
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
	'ToolTipString','',...
	'FontName',PF.times,...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',FS(16),...
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

if strcmp(spm_platform('filesys'),'win')
	off=100;
	uicontrol(F,'Style','PopUp','Tag','DrivesPopup',...
		'ToolTipString','drives',...
		'HorizontalAlignment','Left',...
		'ForegroundColor','r',...
		'String',spm_get('DrivesPullDownStr'),...
		'Callback',['spm_get(''cd''),'...
		    'set(gcbo,''String'',spm_get(''DrivesPullDownStr''))'],...
		'Position',[010+335-off 303 off 022].*WS);
else
	off=0;
end
uicontrol(F,'Style','PopUp','Tag','SubDirsPopup',...
	'ToolTipString','cd to subdirectories of this directory',...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'String','SubDirectories...',...
	'Callback','spm_get(''cd'')',...
	'Position',[010 303 335-off 022].*WS);

uicontrol(F,'Style','PushButton','Tag','CDhome',...
	'String','home',...
	'ToolTipString','cd to your home directory',...
	'ForegroundColor','r',...
	'CallBack',['spm_get(''cd'',''',deblank(LastDirs(2,:)),''')'],...
	'Position',[345 303 045 022].*WS);

uicontrol(F,'Style','Pushbutton',...
	'String','Filter',...
	'ToolTipString','reset filter to ''*''',...
	'CallBack','spm_get(''dir'',[],''*'')',...
	'UserData','*',...
	'Position',[010 276 040 026].*WS);

uicontrol(F,'Style','Pushbutton','String','',...
	'ToolTipString','click to reset filter to initial value',...
	'CallBack',[	'spm_get(''dir'',[],get(findobj(gcbf,''Tag'',',...
			'''Filter''),''UserData''))'],...
	'Position',[048 276 099 026].*WS);

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
	'Position',[150 278 030 022].*WS);

uicontrol(F,'Style','Pushbutton','String','A&D',...
	'ToolTipString','select all items in this directory & Done',...
	'Callback','spm_get(''All'',gcbf);spm_get(''Done'')',...
	'ForegroundColor','m',...
	'Interruptible','off','BusyAction','Cancel',...
	'Position',[180 278 030 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Edit',...
	'ToolTipString','edit list of items',...
	'ForegroundColor','b',...
	'Callback','spm_get(''Edit'')',...
	'Position',[215 278 040 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Keybd',...
	'ToolTipString','switch to keyboard input of items',...
	'ForegroundColor','k',...
	'Callback','spm_get(''GUI2CmdLine'')',...
	'Position',[260 278 040 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Reset',...
	'ToolTipString','reset selection',...
	'ForegroundColor','r',...
	'Callback','spm_get(''Reset'',gcbf)',...
	'Position',[305 278 040 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Done',...
	'ToolTipString','done - press when completed selection of items',...
	'ForegroundColor','m',...
	'Tag','Done','UserData',1,...
	'Callback','spm_get(''Done'')',...
	'Interruptible','off','BusyAction','Cancel',...
	'Position',[350 278 040 022].*WS);

uicontrol(F,'Style','Slider','Tag','ScrollBar',...
	'ToolTipString','scroll',...
	'Callback',[...
		'set(gca,''Units'',''Points''),',...
		'set(gca,''Position'',',...
		'get(gcbo,''UserData'')-[0,get(gcbo,''Value''),0,0])'],...
	'Position',[380 030 020 240].*WS)
	
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
	tmp = max(find(str==filesep));
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
if nargin<5
	Filter=0;
else
	Filter=varargin{5};
end
if nargin<4 Prompt='Select files...'; else Prompt=varargin{4}; end
if nargin<3 n=struct('fi',1,'mn',0,'mx',Inf); else n=varargin{3}; end
if nargin<2 Vis='on'; else Vis=varargin{2}; end

%-Recover SelFileWin figure number
F  = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

cF = get(0,'CurrentFigure');				%-Save current figure

switch lower(Vis), case 'close'
	close(F)
	varargout = {[],cF};
	return
case {'off','reset'}
	varargout = {F,cF};				%-Return figure handles
	if isempty(F), return, end			%-...if no SelFileWin
	varargout = {F,cF};				%-Return figure handles
	set(F,'Visible','off')				%-Make window Invisible
	set(findobj(F,'Tag','Done'),'UserData',1)	%-Set Done UserData to 1
	if strcmp(lower(Vis),'reset')
		set(findobj(F,'Tag','Done'),'UserData',-1)
	end
	delete(get(F,'CurrentAxes'))			%-delete 'dir' axes
	drawnow
	return
case 'on'
	if isempty(F), 	F=spm_get('CreateFig'); end	%-Create SelFileWin
	varargout = {F,cF};				%-Return figure handles
	set(findobj(F,'Tag','Done'),'UserData',0)	%-Init. Done UserData
	delete(get(F,'CurrentAxes'))			%-delete 'dir' axes
	set(findobj(F,'Tag','P'),'UserData','')		%-clear P
	set(findobj(F,'Tag','Prompt'),...
		'String',Prompt,...
		'ToolTipString',Prompt,...
		'UserData',n)
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
	figure(F),drawnow				%-PopUp figure
	set(0,'CurrentFigure',cF)			%-Return to prev. figure

otherwise
	error('Unrecognised ''Vis'' option')
end



case 'statusline'
%=======================================================================
% spm_get('StatusLine',nP,n,F)
if nargin<4, F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
	else, F = varargin{4}; end
if nargin<3, n=get(findobj(F,'Tag','Prompt'),'UserData');else,n=varargin{3}; end
if nargin<2, nP=size(get(findobj(F,'Tag','P'),'UserData'),1);else,nP=varargin{2}; end

str=['Selected ',int2str(nP)];
if (n.mn==n.mx), str=sprintf('%s/%d',str,abs(n.mx)); end
if n.fi
	str = [str,' file'];
	if nP~=1, str=[str,'s']; end
else
	str = [str,' director'];
	if nP==1, str=[str,'y']; else, str=[str,'ies']; end
end

if (n.mn~=n.mx)
	if n.mn==0
		if isfinite(n.mx)
			str=sprintf('%s (max %d)',str,n.mx);
		end
	else
		if isinf(n.mx)
			str=sprintf('%s (min %d)',str,n.mn);
		else
			str=sprintf('%s (min %d, max %d)',str,n.mn,n.mx);
		end
	end
end


if (nP>=n.mn) & (nP<=n.mx)
	str = [str,', press "Done" when finished.'];
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
		set(gcbo,'Value',1)
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
set(findobj(F,'Tag','LastDirsPopup'),...
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

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(F,'Pointer','Watch')
delete(get(F,'CurrentAxes')), drawnow

%-Condition parameters and setup variables
%-----------------------------------------------------------------------
if nargin<4, NoComp=0; else, NoComp=1; end
if nargin<3,
	Filter='';
else,
	Filter=varargin{3};
end
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
Rec = [0.02 0.086 0.92 0.58];
hAxes = axes('Parent',F,'Tag','hAxes',...
	'Position',Rec,...
	'DefaultTextInterpreter','none',...
	'Units','Points',...
	'Visible','off');
y     = get(hAxes,'Position');
y0    = floor(y(3));
dy    = 22;
set(hAxes,'Ylim',[0,y0])


%-List current directory
%-----------------------------------------------------------------------
[Files,Dirs] = spm_list_files_expand(WDir,Filter,ListOpts);
if isempty(Dirs)
	text(0,y0,'Permission denied, or non-existent directory',...
		'Parent',hAxes,'FontWeight','bold','Color','r');
	set(findobj(F,'Tag','ScrollBar'),'Visible','off')
	set(F,'Pointer','Arrow')
	return
else	%-Exclude "." directory from Dirs, exclude ".." if in root
	%-Exclude ".." from DirItems (used when selecting directories)
	Dirs(Dirs(:,1)=='.',:)=[];
	DirItems = strvcat('.',Dirs);
	if ~isroot(WDir), Dirs=strvcat('..',Dirs); end
end

%-Create list of directories in pulldown menu
%-----------------------------------------------------------------------
set(findobj(F,'Tag','SubDirsPopup'),'Value',1,...
	'String',strvcat('SubDirectories...',Dirs))

%-Create list of directories with appropriate 'ButtonDownFcn': -
%-----------------------------------------------------------------------
y     = y0-dy;
if size(Dirs,1)
    for i = 1:size(Dirs,1)
	h = text(0,y,deblank(Dirs(i,:)),...
		'Parent',hAxes,...
		'Tag','DirName',...
		'FontWeight','bold','Color','r',...
		'UserData',deblank(Dirs(i,:)),...
		'ButtonDownFcn','spm_get(''cd'',get(gcbo,''UserData''))',...
		'HandleVisibility','off');
	y = y - dy;
    end
    %-Note lowest position used (in points)
    set(h,'Units','Points')
    my = get(h,'Position')*[0;1;0];
else
    my = 0;
end

%-Files or directories?
%-----------------------------------------------------------------------
if n.fi
	%-Select file(s) - omit ".*" dot files
	if ~isempty(Files), Files(Files(:,1)=='.',:)=[]; end
	Items  = Files;
else	%-Select directory/ies - don't compress the listing
	Items  = DirItems;
	NoComp = 1;
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
	text(0.45,y0,'Summary View - Click to expand',...
		'Parent',hAxes,...
		'Color','w',...
		'FontSize',spm('FontSize',9),...
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
if size(IName,1)
    for i = 1:size(IName,1)
	cIName   = IName(i,IName(i,:)~=' ');
	cItemPos = ItemPos(i,ItemPos(i,:)>0);
	nItems   = length(cItemPos);
	cItems   = deblank(Items(cItemPos,:));
	h = text(0.4,y,cIName,...
		'Parent',hAxes,...
		'Tag','IName',...
		'UserData',cItems,...
		'Color','k',...
		'ButtonDownFcn','spm_get(''Add'')',...
		'HandleVisibility','on');
	if nItems>1;
		h = text(0.39,y-3,int2str(sum(ItemPos(i,:)>0)),...
			'Parent',hAxes,...
			'Color','w',...
			'FontSize',spm('FontSize',8),...
			'FontAngle','Italic',...
			'HorizontalAlignment','right',...
			'UserData',cIName,...
			'ButtonDownFcn',...
			'spm_get(''dir'',[],get(gcbo,''UserData''),1)',...
			'HandleVisibility','off');
	end
	y = y - dy;
    end
    %-Note lowest position used (in points)
    set(h,'Units','Points')
    my = min(my,get(h,'Position')*[0;1;0]);
end

%-Setup scrollbar
%-----------------------------------------------------------------------
h = findobj(F,'Tag','ScrollBar');
if my>0
	set(h,'Visible','off')
else
	set(h,	'UserData',	get(hAxes,'Position'),...
		'Value',	0,...
		'Max',		0,...
		'Min',		my-dy,...
		'SliderStep',	min([0.5,1],[3*dy,(y0-dy)]/(-my+dy)),...
		'Visible',	'on')
end

%-Done
%-----------------------------------------------------------------------
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
fS = char(zeros(nS,Sl));
for i = 1:nS
	c = max(find(S(i,:)~=' '));
	fS(i,:) = S(i,[c:-1:1,Sl:-1:c+1]);
end
varargout = {fS};


case 'filesummary'
%=======================================================================
% [FSpecs,FnamePos]=spm_get('FileSummary',Fnames,Cend,Filter,len)
if nargin<5, len = 3; else len = varargin{5}; end
if nargin<4 | isempty(varargin{4}),
	Filter = '*';
else
	Filter = varargin{4};
end
if nargin<3, Cend='both'; else, Cend = varargin{3}; end
if nargin<2 | isempty(varargin{2}), error('Specify Fnames to summarise')
else Fnames = varargin{2}; end
%-Implicit assumption in this code is that no two files have the same name

if strcmp(Cend,'front')
	if size(Fnames,1)==1 | size(Fnames,2)<len
		%-Only one filename | text shorter than required common length
		varargout = {deblank(Fnames),[1:size(Fnames,1)]'};
		return
	end
	WildCard    = Filter(max(find(Filter=='*')):end);
	nFnames     = size(Fnames,1);
	[sFnames,I] = sortrows(Fnames);
	tmp         = diff(abs(sFnames));
	i           = [0,find(any(tmp(:,1:len),2))',nFnames];
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
		else
			tmp    = min(find([any(diff(abs(Cfnames)),1),1]))-1;
			FSpecs = strvcat(FSpecs,[Cfnames(1,1:tmp),WildCard]);
		end
		%-Note indices of Fnames in this group, retain original ordering
		tI       = sort(I(cI));
		tmp	 = max([nCfnames,size(FnamePos,2)]);
		FnamePos = [FnamePos,...
				zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
				tI', zeros(1,tmp-length(tI))];
	end
	%-Order groups by original ordering of first member occurence
	[null,tmp] = sort(FnamePos(:,1));
	FSpecs     = FSpecs(tmp,:);
	FnamePos   = FnamePos(tmp,:);

elseif strcmp(Cend,'end')
	fFnames           = spm_get('strfliplr',Fnames);
	fFilter           = fliplr(Filter);
	[FSpecs,FnamePos] = spm_get('FileSummary',fFnames,'front',fFilter,len);
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
		    %-Assumming unique input => one blank line at most
		    % NB: This may change ordering of groups by first member
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
		%-Look for common endings
		[tSpecs,tI] = spm_get('FileSummary',Cfnames,'end',Filter,len);
		%-If no common endings then don't expand to individual files!
		if size(tI,1)>1 & size(tI,2)==1
			tSpecs = '*';
			tI = tI';
		end
		for j = 1:size(tSpecs,1)
		    cI       = chI(tI(j,tI(j,:)>0));
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
	if (nP+nItems>n.mx), break, end

	%-Add items to P
	if WDir(end)~=filesep, WDir = [WDir,filesep]; end
	FPath  = [repmat(WDir,nItems,1),Items];
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

%-Return new P to holding object
%-----------------------------------------------------------------------
set(findobj(F,'Tag','P'),'UserData',deblank(P));

%-Update StatusLine
%-----------------------------------------------------------------------
spm_get('StatusLine',size(P,1),n,F)


case 'delete'
%=======================================================================
% spm_get('Delete',h)
F      = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
if nargin<2, h=gcbo; else, h=varargin{2}; end
WDir   = get(findobj(F,'Tag','WDir'),'UserData');
if WDir(end)~=filesep, WDir = [WDir,filesep]; end
P      = get(findobj(F,'Tag','P'),'UserData');
n      = get(findobj(F,'Tag','Prompt'),'UserData');
IName  = get(h,'String');
IName(1:find(IName==':')) = [];
Items  = get(h,'UserData');
nItems = size(Items,1);
FPath  = [repmat(WDir,nItems,1),Items];
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
if nargin<2, n=struct('fi',1,'mn',0,'mx',Inf); else n=varargin{2}; end

if n.mx==0, varargout={''}; return, end
nP=size(P,1);

%-Print prompt
%-----------------------------------------------------------------------
fprintf('\n'), fprintf('%c',repmat('=',1,70)), fprintf('\n')
fprintf('\t%s\n',Prompt)
fprintf('%c',repmat('=',1,70)), fprintf('\n')
fprintf('Current working directory: %s\n',WDir)
fprintf('(Prepended to relative paths)\n')
if nP>0
	fprintf('\nAlready selected:\n')
	for pP = [1:size(P,1)], fprintf('\t%s\n',P(pP,:)), end
end

if nP>=n.mx
	fprintf('\nMax # items already selected: Returning...\n\n' )
	varargout={P};
	return
end

%-Print prompts for input
%-----------------------------------------------------------------------
fprintf('\nEnter paths ')
if n.fi
	str = 'file'; if n.mx>1, str=[str,'s']; end
else
	str = 'director'; if n.mx==1, str=[str,'y']; else, str=[str,'ies']; end
end
if AllowEnd | n.mn==0
	if isfinite(n.mx)
		fprintf('to at most %d %s:\n',n.mx,str)
	else
		fprintf('to %s:\n',str)
	end
else
	if (n.mn==n.mx)
		fprintf('to %d %s:\n',n.mx,str)
	else
		if isinf(n.mx)
			fprintf('to at least %d %s:\n',n.mn,str)
		else
			fprintf('to %d-%d %s:\n',n.mn,n.mx,str)
		end
	end
end

%-AllowEnd if asked for, or if a range of item numbers is OK
%-----------------------------------------------------------------------
Tstr = 'END';
if AllowEnd | (n.mn~=n.mx)
	fprintf(' (Type "%s" to terminate input.)\n',Tstr)
end

%-Get files
%-----------------------------------------------------------------------
Done=0;
while ~( Done | (nP>=n.mx) )
	%-Get input
	str = []; while isempty(str)
		str=ddeblank(input(sprintf('  %3d  : ',nP+1),'s'));end
	%-Prepend WDir to incomplete pathnames
	if (~isabspath(str))&(~strcmp(str,Tstr)), str=fullfile(WDir,str); end
	if n.fi,OK=exist(str)==2;else,OK=exist(str)==7;end
	OKend = AllowEnd | (nP>=n.mn&nP<=n.mx);
	while ~( OK | (strcmp(str,Tstr) & OKend) )
		if strcmp(str,Tstr)
			fprintf('%c\tSelect %d-%d items!\n',7,n.mn,n.mx)
			else, fprintf('%c\t%s isn''t valid!\n',7,str), end
		str=[]; while isempty(str)
			str=ddeblank(input(sprintf('  %3d  : ',nP+1),'s'));end
		if (~isabspath(str))&(~strcmp(str,Tstr))
			str = fullfile(WDir,str); end
		if n.fi,OK=exist(str)==2;else,OK=exist(str)==7;end
		OKend = AllowEnd | (nP>=n.mn&nP<=n.mx);
	end
	if ~strcmp(str,Tstr), P=strvcat(P,str); end
	nP = nP+1;
	Done = strcmp(str,Tstr);
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

h = uicontrol(F,'Style','Text',...
	'String',Prompt,...
	'ToolTipString',Prompt,...
	'FontName',spm_platform('font','times'),...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',spm('FontSize',16),...
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
	'FontName',spm_platform('font','times'),...
	'FontWeight','Bold',...
	'FontAngle','Italic',...
	'FontSize',spm('FontSize',16),...
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
	'String',cellstr(char(get(findobj(F,'Tag','P'),'UserData'))),...
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

str=[]; if isfinite(n.mx), str=['Use only top ',int2str(n.mx),' lines.  ']; end
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
	P = P(~strcmp(P,''));
	if (isfinite(n.mx))&(~isempty(P)), P=P(1:min(length(P),n.mx)); end
	%-Write P to 'P' 'Tag'ged object
	set(findobj(F,'Tag','P'),'UserData',char(P));
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
nP = size(get(findobj(F,'Tag','P'),'UserData'),1);
if (nP<n.mn | nP>n.mx)
	if (n.mn==n.mx)
		str = sprintf('Must select %d item(s)!',n.mx);
	else
		str = sprintf('Must select %d-%d items!',n.mn,n.mx);
	end
	spm('Beep')
	msgbox(str,sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'error','non-modal')
	drawnow

	% Not sure about this bit - JA
	% set(F,'windowstyle','modal')

	return
end

%-Done, set Done UserData tag for handling
set(findobj(F,'Tag','Done'),'UserData',1)


case 'cpath'
%=======================================================================
% cpath = spm_get('cpath',path,cwd)
if nargin<3, cwd=pwd; else, cwd=spm_get('CPath',varargin{3}); end
if nargin<2, error('insufficient arguments'), end
cpath = cellstr(varargin{2});

sepchar = filesep;			%-File seperator character
dubsep  = [sepchar,sepchar];		%-Double file separator character
rootlen = spm_platform('rootlen');	%-Length of root paths ('/' & 'C:\' &c)

for i=1:prod(size(cpath))

	%-Prepend cwd to relative pathnames
	%---------------------------------------------------------------
	%-Try to figure out "~username" paths on unix/linux
	%if ~isabspath(cpath{i}), cpath{i}=fullfile(cwd,cpath{i}); end
	if ~isabspath(cpath{i}),
		if isunix & ~isempty(cpath{i}) & cpath{i}(1)=='~'
			try
				opwd=pwd;
				cd(cpath{i});
				cpath{i}=pwd;
				cd(opwd);
			catch
				cpath{i}=fullfile(cwd,cpath{i});
			end;
		else
			cpath{i}=fullfile(cwd,cpath{i});
		end;
	end

	%-Sort out stationary relative pathnames './' & '/.'
	%---------------------------------------------------------------
	%-Remove midpath '/./'
	cpath{i}=strrep(cpath{i},[sepchar,'.',sepchar], sepchar);
	
	%-Remove '.' from trailing '/.'
	if length(cpath{i})>=2 & strcmp(cpath{i}(end-1:end),[sepchar,'.'])
		cpath{i}(end)=''; end
	
	%-Remove trailing '/'s (unless root, in which case ensure trailing '/')
	if isroot(cpath{i})
		if cpath{i}(end)~=sepchar, cpath{i}(end+1)=sepchar; end
	else
		while(length(cpath{i})>rootlen) & cpath{i}(end)==sepchar
			cpath{i}(end)='';
		end
	end
	
	%-Remove any residual '//'s
	while length(cpath{i})>=2 & length(findstr(cpath{i},dubsep))
		cpath{i}=strrep(cpath{i},dubsep,sepchar); end

	%-Sort out relative pathnames '/..'
	%---------------------------------------------------------------
	%-Process midpath '/../'
	t0 = 0;
	downstr = [sepchar,'..',sepchar];
	while length(cpath{i})>=rootlen+4 & max([0,findstr(cpath{i},downstr)])>t0
		t2 = t0+min(findstr(cpath{i}(t0+1:end),downstr));
		t1 = t0+max([0,find(cpath{i}(t0+1:t2-1)==sepchar)]);
		if t1==t0 | strcmp(cpath{i}(t1:t2),downstr)
			t0=t2;
		else
			cpath{i}(t1:t2+2)='';
		end
	end

	%-Process trailing '/..' (Don't remove first '/')
	if length(cpath{i})>rootlen+2 & ...
	   strcmp(cpath{i}(end-2:end),[sepchar,'..'])
		t2 = length(cpath{i})-2;
		t1 = t0+max([0,find(cpath{i}(t0+1:t2-1)==sepchar)]);
		if t1~=t0 & ~strcmp(cpath{i}(t1:t2),downstr)
			cpath{i}(max(t1,rootlen+1):t2+2)='';
		end
	end

	%-Remove any leading '/../'s
	while length(cpath{i})>=rootlen+3 & ...
	      strcmp(cpath{i}(rootlen+[0:3]),downstr)
		cpath{i}(rootlen+[1:3])='';
	end

	%-Canonicalise '/..' to '/'
	if length(cpath{i})==rootlen+2 & ...
	   strcmp(cpath{i}(rootlen+[0:2]),[sepchar,'..'])
		cpath{i}(rootlen+[1:2])='';
	end

	%-Remove any residual '//'s
	while length(cpath{i})>=2 & length(findstr(cpath{i},dubsep))
		cpath{i}=strrep(cpath{i},dubsep,sepchar); end


	%-Final canonicalisations for relative pathnames
	%---------------------------------------------------------------
	% %-Remove leading './'
	% if length(cpath{i})>2 & strcmp(cpath{i}(1:2),['.',sepchar])
	% 	cpath{i}(1:2)=''; end
	
	%-Canonicalise leading './../' to '../'
	if length(cpath{i})>5 & strcmp(cpath{i}(1:5),['.',downstr])
		cpath{i}(1:2)=''; end
	%-Canonicalise './..' to '..'
	if strcmp(cpath{i},['.',sepchar,'..']), cpath{i}='..'; end

end

if ischar(varargin{2}), varargout={char(cpath)}; else, varargout={cpath}; end


case 'drivespulldownstr'
%=======================================================================
% str = spm_get('DrivesPullDownStr')
drivestr  = spm_win32utils('drives');
n         = length(drivestr);
varargout = ...
	{['Drives...',reshape(['|'*ones(1,n);drivestr;':'*ones(1,n)],1,3*n)]};


case 'files'
%=======================================================================
% FORMAT [P,dir] = spm_get('Files',dir,fil)
if nargin<3,
	fil='*';
else,
	fil = varargin{3};
end
if nargin<2, dir='.'; else, dir = varargin{2}; end
dir = spm_get('CPath',dir);
[P,null] = spm_list_files_expand(dir,fil,ListOpts);

%-Make into full pathnames if nargout<2
if ~isempty(P) & nargout<2
	nP = size(P,1);
	P = [repmat(dir,nP,1), repmat(filesep,nP,1), P];
end

varargout = {P,dir};


otherwise
%=======================================================================
error('Illegal Action string')


%=======================================================================
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function so = ddeblank(si)
%=======================================================================
if nargin==0, error('insufficient arguments'), end
if isempty(si), so=''; return, end
if ~ischar(si), error('input must be a string'), end

b = any((si~=' ' & si~=0),1);
if all(~b), so=''; return, end
so = si(:,max(1,min(find(b))):min(max(find(b)),length(b)));


function rootf = isroot(path)             %%-Platform specifics (MBrett)
%=======================================================================
%-Returns true if path is root directory, else false
switch (spm_platform('filesys'))
case 'unx'
	rootf = strcmp(path,'/');
case 'win'
	%-Accepts e.g e: and e:\ as root paths
	if (length(path)==2 & path(2)==':') | ...
	   (length(path)>2 & strcmp(path(2:end),':\'))
		rootf = 1;
	else
		rootf = 0;
	end
otherwise
	error('isroot not coded for this filesystem');
end


function absf = isabspath(path)           %%-Platform specifics (MBrett)
%=======================================================================
%-Returns true if path is absolute, false if relative (or empty)
switch (spm_platform('filesys'))
case 'unx'
	if (~isempty(path) & path(1)=='/'), absf=1; else, absf=0; end
case 'win'
	if (length(path)>1 & path(2)==':'), absf=1; else, absf=0; end
otherwise
	error('isabspath not coded for this filesystem');
end
return

