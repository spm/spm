function [R1,R2] = spm_get(Action,P2,P3,P4,P5,P6)
% User interface : filename selection
% FORMAT [P] = spm_get(n,Suffix,Prompt,Prefix,NewWDir,CmdLine);
% n          - 0 - returns with empty P
%            - Positive integer - prompts for n filepaths
%            - +Inf - prompts for filepaths until "Done" pressed
%            - Negative integer - prompts for n directories
%            - -Inf - prompts for directories until "Done" pressed
% Suffix     - Filename filter {i.e. *Suffix} 
%            - numeric 0 for previous filtering
% Prompt     - Prompt string
% Prefix     - Prefix Filename filter {i.e. Prefix*}
% NewWDir    - New working directory
% CmdLine    - CmdLine override switch, 0 for GUI, 1 for Command line.
%
% P          - string matrix of full pathnames {one string per row}
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
% The filter string for filenames is built as Prefix*Suffix, and
% appears in a single editable window. Subsequent calls to spm_get
% callbacks with the same filter string do not change the value of the
% filter string, allowing the users edits of the filer string to
% persist. Clicking the Filter button resets the filter to the original
% filter string passed by the program. Clicking the bar at the
% right-hand end of the filter window sets the filter string to '*',
% resulting in all files/directories being displayed, in summary view
% if there are sufficient files. (Bars above and below the filter
% window (if available) implement a simple "memory" for filter
% strings.  The top bar stores the current filter string, the lower bar
% recalls it.)
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
% FORMAT R1 = spm_get(Action,P2,P3,P4,P5,P6)
%           - A multi function function!
%
% FORMAT F = spm_get('CreateFig',LastDirs,Filter)
% Sets up the file selection window as next avaliable figure.
% F        - Figure used.
% LastDirs - String matrix for directory history.
% Filter   - Filename filter.
%
% FORMAT spm_help('FigKeyPressFcn',ch)
% Callback for handling keyboard accelerators
%
% FORMAT F = spm_get('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.
% F        - Figure used.
% Vis      - Figure visibility (string, 'on' or 'off')
% n        - Number of filepaths/directories to obtain
% Prompt   - Prompt string.
% Filter   - Filename filter.
% WDir     - [optional] New Working directory
%
% FORMAT spm_get('StatusLine',nP,n)
% Updates FileSelWin status line
% nP       - #Items already selected
% n        - #Items to select
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
% Filter   - Filename filter - Suffix filter is appended to 'front' items
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
% FORMAT spm_get('Reset')
% Resets the selection and redisplays the directory
%
% FORMAT spm_get('All')
% Adds all remaining unselected items to the selection
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
%_______________________________________________________________________
% Andrew Holmes


%-Parameters
%=======================================================================
%-Default to unlimited file get when no arguments
if nargin<1 Action=+Inf; end
%-Allow pointer jumping?
PJump = 0;

if ~isstr(Action)
%=======================================================================
% P = spm_get(n,Suffix,Prompt,Prefix,NewWDir,CmdLine)
%fprintf('Doing: spm_get(...\n')

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<6 CmdLine=[]; else CmdLine=P6; end
if nargin<5 NewWDir=''; else NewWDir=P5; end
if nargin<4 Prefix=''; else Prefix=P4; end
if nargin<3 Prompt='Select files...'; else Prompt=P3; end
if nargin<2 Suffix=[]; else Suffix=P2; end
if isempty(Suffix), Suffix=0; end
n = Action;

if (n==0), P=[]; return, end

if isempty(CmdLine)
	global CMDLINE
	if ~isempty(CMDLINE), CmdLine = CMDLINE; else, CmdLine=0; end
end

%-NB: spm_get callbacks use single filter string, or 0 for previous filtering
if isstr(Suffix)
	if Suffix(1)~='*', Suffix = ['*',Suffix]; end
	Filter = [Prefix,Suffix];
else
	Filter=0;
end

%-Computation
%=======================================================================

if CmdLine

	%-Command line version, if global CMDLINE is true, & not overridden
	%---------------------------------------------------------------
	P=spm_get('CmdLine',n,Prompt,[],NewWDir);

else

	%-GUI version
	%---------------------------------------------------------------
	%-Set up FileSelWin & JumpCursor
	[F,cF] = spm_get('Initialise','on',n,Prompt,Filter,NewWDir);
	PLoc = get(0,'PointerLocation');
	FRec = get(F,'Position');
	if PJump, set(0,'PointerLocation',...
		[FRec(1)+FRec(3)/2, FRec(2)+FRec(2)/2]), end
	
	%-Wait until filenames have been selected
	hDone = findobj(F,'Tag','Done');
	waitfor(hDone,'UserData')
	
	%-Recover P
	P=get(findobj(F,'Tag','P'),'UserData');
	
	%-Reset and hide SelFileWin
	spm_get('Initialise','off');

	%-Return focus to previous figure (if any), & put cursor back
	set(0,'CurrentFigure',cF)
	%-Jump cursor back to previous location
	if PJump, set(0,'PointerLocation',PLoc), end

end % (if CmdLine)

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_get : ',Prompt,' :'],P); end

R1 = P;
return

elseif strcmp(Action,'CreateFig')
%=======================================================================
% F = spm_get('CreateFig',LastDirs,Filter)
%fprintf('Doing: spm_get(''CreateFig''...\n')

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<3 Filter=[]; else Filter=P3; end
if nargin<2 LastDirs=[]; else LastDirs=P2; end
if isempty(LastDirs), LastDirs=pwd; end
LastDirs=str2mat(LastDirs,getenv('HOME'));
if (exist('spm.m')==2), LastDirs=str2mat(LastDirs,spm('Dir')); end
	

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
	'DefaultTextFontSize',2*round(12*min(WS)/2),...
	'DefaultUicontrolFontSize',2*round(12*min(spm('GetWinScale'))/2),...
	'DefaultUicontrolInterruptible','on',...
	'Visible','off');
R1=F;

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
	'String',str2mat('Previous Directories...',LastDirs),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'UserData',size(LastDirs,1),...
	'Callback','spm_get(''cd'')',...
	'Position',[010 347 335 022].*WS);

uicontrol(F,'Style','PushButton','Tag','CDpwd',...
	'String','pwd',...
	'ForegroundColor','r',...
	'UserData',0,...
	'CallBack','spm_get(''cd'',pwd)',...
	'Position',[345 347 045 022].*WS);

WDir=deblank(LastDirs(1,:));
uicontrol(F,'Style','Edit','String',WDir,...
	'Tag','WDir','UserData',WDir,...
	'ForegroundColor','r','BackgroundColor',[.8,.8,1],...
	'HorizontalAlignment','Left',...
	'CallBack','spm_get(''cd'')',...
	'Position',[010 325 380 022].*WS);

uicontrol(F,'Style','PopUp','Tag','SubDirsPopup',...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'String','SubDirectories...',...
	'Callback','spm_get(''cd'')',...
	'Position',[010 303 335 022].*WS);

uicontrol(F,'Style','PushButton','Tag','CDhome',...
	'String','home',...
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
	'CallBack',[...
		'spm_get(''dir'',[],get(findobj(gcbf,''Tag'',',...
			'''Filter''),''UserData''))'],...
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
	'CallBack','spm_get(''dir'',[],''*'')',...
	'Position',[145 278 005 022].*WS);

uicontrol(F,'Style','Edit','Tag','Filter',...
	'String','*',...
	'UserData','*',...
	'BackgroundColor',[.8,.8,1],...
	'Callback','spm_get(''dir'')',...
	'Position',[050 278 095 022].*WS);

uicontrol(F,'Style','Pushbutton','String','All',...
	'Callback','spm_get(''All'')',...
	'Position',[155 278 030 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Edit',...
	'ForegroundColor','b',...
	'Callback','spm_get(''Edit'')',...
	'Position',[190 278 040 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Keybd',...
	'ForegroundColor','k',...
	'Callback','spm_get(''GUI2CmdLine'')',...
	'Position',[235 278 050 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Reset',...
	'ForegroundColor','r',...
	'Callback','spm_get(''Reset'')',...
	'Position',[290 278 050 022].*WS);

uicontrol(F,'Style','Pushbutton','String','Done',...
	'ForegroundColor','m',...
	'Tag','Done','UserData',0,...
	'Callback','spm_get(''Done'')',...
	'Interruptible','off','BusyAction','Cancel',...
	'Position',[345 278 045 022].*WS);

uicontrol(F,'Style','Pushbutton','String','',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''UserData''))'],...
	'Position',[370 265 020 005].*WS);

uicontrol(F,'Style','Pushbutton','String','/\',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''Position'')-[0 0.5 0 0])'],...
	'Position',[370 245 020 020].*WS);

uicontrol(F,'Style','Pushbutton','String','\/',...
	'Callback',[...
		'set(gca,''Units'',''Normalized'');',...
		'set(gca,''Position'',get(gca,''Position'')+[0 0.5 0 0])'],...
	'Position',[370 040 020 020].*WS);

uicontrol(F,'Style','Frame','Tag','StatusArea',...
	'Position',[001 000 400 030].*WS);

if exist('spm_help.m')==2
	uicontrol(F,'Style','Pushbutton','String','?',...
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
	'spm_get(''FigKeyPressFcn'',get(gcf,''CurrentCharacter''))')
return


elseif strcmp(Action,'FigKeyPressFcn')
%=======================================================================
% spm_help('FigKeyPressFcn',ch)
if nargin<2, ch=get(gcbf,'CurrentCharacter'); else, ch=P2; end

%-Empty ch - shift/control/&c.
if isempty(ch), return, end

%-Keyboard accelerators for item selection
%-----------------------------------------------------------------------
if any(abs(ch)==[10,11,19]) % ^J,^K,^S
	H = flipud(findobj(get(gca,'Children'),...
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
hCDpwd    = findobj(gcbf,'Tag','CDpwd');
EditWDir  = get(hCDpwd,'UserData');
if EditWDir, h = findobj(gcbf,'Tag','WDir');
	else, h = findobj(gcbf,'Tag','Filter'); end
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
return


elseif strcmp(Action,'Initialise')
%=======================================================================
% [F,cF] = spm_get('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.
%fprintf('Doing: spm_get(''Initialise''...\n')

if nargin<6 WDir=''; else WDir=P6; end
if nargin<5 Filter=0; else Filter=P5; end
if nargin<4 Prompt='Select files...'; else Prompt=P4; end
if nargin<3 n=Inf; else n=P3; end
if nargin<2 Vis='on'; else Vis=P2; end

%-Save current figure
cF = get(0,'CurrentFigure');

%-Recover SelFileWin figure number
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

%-If closing SelFileWin then close and return
if strcmp('Vis','Close'), if ~isempty(F), close(F), end, R1=[]; return, end

%-If no SelFileWin, then create one
if isempty(F), F=spm_get('CreateFig'); end

%-delete 'dir' text axes
delete(gca)

%-If Vis=='off', make invisible
if strcmp(Vis,'off'), set(F,'Visible','off'), end

%-clear P
set(findobj(F,'Tag','P'),'UserData','')

%-Reset prompt, & n, stored in UserData
set(findobj(F,'Tag','Prompt'),'String',Prompt,'UserData',n)

%-Reset StatusLine
spm_get('StatusLine',0,n)

%-Done UserData to zero
set(findobj(F,'Tag','Done'),'UserData',0)

%-KeyPressFcn to edit WDir widget initially
set(findobj(F,'Tag','CDpwd'),'UserData',1)

%-Filter string
%-Only set if different to previously passed one, held as UserData of
% Filter Edit uicontrol.
if isstr(Filter)
	if ~strcmp(Filter,get(findobj(F,'Tag','Filter'),'UserData'))
	    set(findobj(F,'Tag','Filter'),'String',Filter,'UserData',Filter)
	end
end

%-Change to new working directory if one was specified
if ~isempty(WDir), spm_get('cd',WDir,1), end

%-Make current & visible and do 'dir', if Vis='on'
if strcmp(Vis,'on')
	spm_get('dir')
	figure(F)
	set(F,'Visible','on')
end

%-Return figure handles
R1 = F; R2 = cF;

return


elseif strcmp(Action,'StatusLine')
%=======================================================================
% spm_get('StatusLine',nP,n)
%fprintf('Doing: spm_get(''StatusLine''...\n')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
if nargin<3, n=get(findobj(F,'Tag','Prompt'),'UserData'); else, n=P3; end
if nargin<2, nP=size(get(findobj(F,'Tag','P'),'UserData'),1); else, nP=P2; end

str=['Selected ',int2str(nP)];
if finite(n), str=[str,'/',int2str(abs(n))]; end
if (n>=0)
	str = [str,' file'];
	if nP~=1, str=[str,'s']; end
else
	str = [str,' director'];
	if nP==1, str=[str,'y']; else, str=[str,'ies']; end
end % (if)

if isinf(n)
	str = [str,', press "Done" when finished.'];
elseif nP==abs(n)
	str = [str,', press "Done".'];
else
	str = [str,'.'];
end

set(findobj(F,'Tag','StatusLine'),'String',str)
return


elseif strcmp(Action,'cd')
%=======================================================================
% spm_get('cd',NewDir,bDirList)
%fprintf('Doing: spm_get(''cd''...\n')

%-Condition arguments
%-----------------------------------------------------------------------
if (nargin<3), bDirList=1; else, bDirList = ~P3; end
if (nargin<2), NewDir=''; else, NewDir=P2; end

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
	NewDir=P2;
end % (if)
if isempty(NewDir), NewDir=WDir; end
NewDir=deblank(NewDir);

%-Condition relative pathnames
%-----------------------------------------------------------------------
if NewDir(1)~='/'
	if NewDir(1)~='.'
		%-Prepend current working directory
		if strcmp(WDir,'/') NewDir=['/',NewDir];
			else NewDir=[WDir,'/',NewDir]; end
	else
		%-Sort out "dot" options...
		NewDir(1)=[];
		if ~isempty(NewDir) & NewDir(1)=='.'
			%-NewDir started "..", go down a directory from WDir
			WDir=WDir(1:max(find(WDir=='/'))-1);
			if isempty(WDir), WDir='/'; end
			NewDir(1)=[];
		end
		NewDir=[WDir,NewDir];
	end % (if NewDir(1)~='.')
end % (if NewDir...)

if length(NewDir)>=3
	%-Sort out trailing '/..'
	tmp = fliplr(NewDir);
	if all(tmp(1:3)=='../')
		tmp(1:3)='';
		if ~isempty(tmp), tmp(1:min(find(tmp=='/')))=''; end
		if ~isempty(tmp)
			NewDir=fliplr(tmp);
		else
			NewDir='/';
		end
	end
end

% if any(findstr('..',NewDir)
% 	%-Sort out any '..'s in the directory path
% end

if length(NewDir)>=2
	%-Sort out trailing '.' & '/'
	tmp = fliplr(NewDir);
	if all(tmp(1:2)=='./')
		tmp(1:2)='';
		if ~isempty(tmp)
			NewDir=fliplr(tmp);
		else
			NewDir='/';
		end
	elseif tmp(1)=='/'
		NewDir=NewDir(1:length(NewDir)-1);
	end
end


%-Changing directory, delete current directory listing (Done again in
% Action=='dir', but neater to do here in advance on this occasion)
%-----------------------------------------------------------------------
set(0,'CurrentFigure',F), delete(gca)


%-Set up LastDirs
%-----------------------------------------------------------------------
%-Recover LastDirs & size of fixed part from popup menu object
LastDirs=get(findobj(F,'Tag','LastDirsPopup'),'String'); LastDirs(1,:)=[];
Fs=get(findobj(F,'Tag','LastDirsPopup'),'UserData');

%-Add new directory
LastDirs=str2mat(NewDir,LastDirs);

%-Extract fixed and variable parts of LastDirs
FLastDirs=LastDirs(size(LastDirs,1)-[Fs-1:-1:0],:);
LastDirs(size(LastDirs,1)-[Fs-1:-1:0],:)=[];

%-Delete any replications of NewDir within the variable part of LastDirs
if size(LastDirs,1)>1
	IDRows=[0, all( LastDirs(2:size(LastDirs,1),:)'==...
		setstr(ones(size(LastDirs,1)-1,1)*LastDirs(1,:))' ) ];
	LastDirs(find(IDRows),:)=[];
end % (if)

%-Delete NewDir from top of LastDirs if present in fixed part of LastDirs
if any( all(FLastDirs'==setstr(ones(size(FLastDirs,1),1)*LastDirs(1,:))') )
	LastDirs(1,:)=[]; end

%-Lose directories from the bottom of variable LastDirs if oversize
Vs = 10; %-Maximum allowable number of variable LastDirs entries
if size(LastDirs,1)>Vs LastDirs(Vs+1:size(LastDirs,1),:)=[]; end

%-Put LastDirs and FLastDirs back together
LastDirs = [LastDirs; FLastDirs];


%-Write LastDirs to LastDirsPopup, write new directory to WDir Edit widget
%-----------------------------------------------------------------------
set(findobj(F,'Tag','LastDirsPopup'),'Value',1,...
	'String',str2mat('Previous Directories...',LastDirs))
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
return


elseif strcmp(Action,'dir')
%=======================================================================
% spm_get('dir',WDir,Filter,NoComp)
%fprintf('Doing: spm_get(''dir''...\n')
%
% Creates list of text objects with an associated
% ButtonDownFcn functions that call spm_get callbacks for processing
%
%-WDir - Directory, defaults to UserData of 'WDir' tagged object in
% SelFileWin figure

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(0,'CurrentFigure',F)
set(F,'Pointer','Watch')
delete(gca), drawnow

%-Condition parameters and setup variables
%-----------------------------------------------------------------------
if nargin<4, NoComp=0; else, NoComp=1; end
if nargin<3, Filter=''; else, Filter=P3; end
if nargin<2, WDir=''; else, WDir = P2; end
if isempty(WDir), WDir = get(findobj(F,'Tag','WDir'),'UserData');
	else, WDir = deblank(WDir); end

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
hAxes = axes('Position',Rec,...
	'DefaultTextInterpreter','none',...
	'DefaultTextHandleVisibility','off',...
	'Units','Points',...
	'Visible','off',...
	'UserData',Rec,...
	'Tag','hAxes');
y     = floor(get(hAxes,'Position'));
y0    = y(3);
dy    = 22;
set(hAxes(1),'Ylim',[0,y0])


%-List current directory
%-----------------------------------------------------------------------
%fprintf('Doing spm_get(''dir'':\tspm_list_files...\n')
[Files,Dirs] = spm_list_files(WDir,Filter);
if isempty(Dirs)
	text(0,y0,'Permission denied, or non-existent directory',...
		'FontWeight','bold','Color','r');
end
%-Lose dot files & directories, only show .. as a dir if not in /
if ~isempty(Files), Files(Files(:,1)=='.',:)=[]; end
Dirs(Dirs(:,1)=='.',:)=[];
if ~strcmp(WDir,'/'), Dirs=strvcat('..',Dirs); end

%-Create list of directories in pulldown menu
%-----------------------------------------------------------------------
set(findobj(F,'Tag','SubDirsPopup'),'Value',1,...
	'String',str2mat('SubDirectories...',Dirs))

%-Create list of directories with appropriate 'ButtonDownFcn': -
% NB: Gave MatLab Bus errors (ML5.0.0.4064, Solaris2.5.1)
%-----------------------------------------------------------------------
%fprintf('Doing: spm_get(''dir'':\tdirectory list\n')
y     = y0-dy;
for i = 1:size(Dirs,1)
	text(0,y,deblank(Dirs(i,:)),'Tag','DirName',...
		'FontWeight','normal','Color','r',...
		'UserData',deblank(Dirs(i,:)),...
		'ButtonDownFcn','spm_get(''cd'',get(gcbo,''UserData''))',...
		'HandleVisibility','off');
	y = y - dy;
end

%-Files or directories (n<0)
%-----------------------------------------------------------------------
Items = Files; if (n<0) Items = Dirs; end

%-Compressed summary view, or full view?
%-----------------------------------------------------------------------
%fprintf('Doing: spm_get(''dir'':\tcompressed / summary view?\n')
if ~NoComp & size(Items,1)>48
	%-Use a compressed summary view
	if Filter(end)=='*'
		[IName,ItemPos] = spm_get('FileSummary',Items);
	else
		[IName,ItemPos] = spm_get('FileSummary',Items,'front',Filter);
	end
	text(0.5,y0,'Summary View - Click to expand',...
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
%fprintf('Doing: spm_get(''dir'':\titem list\n')
y     = y0-dy;
for i = 1:size(IName,1)
	cIName   = IName(i,IName(i,:)~=' ');
	cItemPos = ItemPos(i,ItemPos(i,:)>0);
	nItems   = length(cItemPos);
	cItems   = Items(cItemPos,:);
	%-Next line strips off redundant space at end of string matrix
	% (The str2mat bit ensures single strings are handled by column)
	cItems(:,all(str2mat(cItems,' ')==' '))=[];
	text(0.35,y,cIName,...
		'Tag','IName',...
		'UserData',cItems,...
		'Color','k',...
		'ButtonDownFcn','spm_get(''Add'')',...
		'HandleVisibility','on');
	if nItems>1;
		text(0.34,y-3,int2str(sum(ItemPos(i,:)>0)),...
			'Color','w',...
			'FontSize',9,...
			'FontAngle','Italic',...
			'HorizontalAlignment','right',...
			'UserData',cIName,...
			'ButtonDownFcn',...
			'spm_get(''dir'',[],get(gco,''UserData''),1)',...
			'HandleVisibility','off');
	end
	y = y - dy;
end
set(F,'Pointer','Arrow')
return


elseif strcmp(Action,'StrSort')
%=======================================================================
% [sS,I]=spm_get('StrSort',S)
%fprintf('Doing: spm_get(''StrSort''...\n')
% Utility string sorting routine for string matrices
if nargin<2, error('Sort what?'), else, S = P2; end

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
R1 = sS;
R2 = I;
return


elseif strcmp(Action,'strfliplr')
%=======================================================================
% fS = spm_get('strfliplr',S)
%fprintf('Doing: spm_get(''strfliplr''...\n')
if nargin<2, error('Flip what?'), else, S = P2; end

%-FlipLR the Fnames matrix, but watch out for trailing blanks/spaces
[nS,Sl] = size(S);
fS = ones(nS,Sl);
for i = 1:nS
	c = max(find(S(i,:)~=' '));
	fS(i,:) = S(i,[c:-1:1,Sl:-1:c+1]);
end
R1 = fS;
return


elseif strcmp(Action,'FileSummary')
%=======================================================================
% [FSpecs,FnamePos]=spm_get('FileSummary',Fnames,Cend,Filter,len)
%fprintf('Doing: spm_get(''FileSummary''...\n')
if nargin<5, len = 3; else len = P5; end
if nargin<4, Filter = ''; else Filter = P4; end
if isempty(Filter), Filter='*'; end
if nargin<3, Cend='both'; else, Cend = P3; end
if nargin<2, error('Specify Fnames to summarise'), else Fnames = P2; end
%-Implicit assumption in this code is that no two files have the same name

if isempty(Fnames), error('Null Fnames'), end

if strcmp(Cend,'front')
	if size(Fnames,1) == 1
		%-Only one filename!
		R1 = deblank(Fnames);
		R2 = 1;
		return
	end
	WildCard    = Filter(max(find(Filter=='*')):length(Filter));
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
			FSpecs = str2mat(FSpecs,Cfnames(1,1:tmp));
		elseif nCfnames == 2
			tmp    = min(find([diff(abs(Cfnames)),1]))-1;
			FSpecs = str2mat(FSpecs,[Cfnames(1,1:tmp),WildCard]);
		else
			tmp    = min(find([any(diff(abs(Cfnames))),1]))-1;
			FSpecs = str2mat(FSpecs,[Cfnames(1,1:tmp),WildCard]);
		end
		tI       = I(cI);
		tmp	 = max([nCfnames,size(FnamePos,2)]);
		FnamePos = [FnamePos,...
				zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
				tI', zeros(1,tmp-length(tI))];
	end
	FSpecs(1,:)   = [];
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
		FSpecs   = str2mat(FSpecs,hSpecs(i,hSpecs(i,:)~=' '));
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
		% Append blank line to prevent any working on a vector
		bI = find(all(str2mat(Cfnames',' ')==' '));
		if ~isempty(bI)
		    % There'll be one blank line at most
		    FSpecs   = str2mat(FSpecs,...
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
		    FSpecs   = str2mat(FSpecs,cFSpec);
		    tmp      = max([length(cI),size(FnamePos,2)]);
		    FnamePos = [FnamePos,...
			zeros(size(FnamePos).*[1 -1]+[0,tmp]);...
			cI, zeros(1,tmp-length(cI))];
		end
	    end
	end
	FSpecs(1,:)   = [];
else
	error('Invalid Cend specifier')
end
R1 = FSpecs;
R2 = FnamePos;
return


elseif strcmp(Action,'Add')
%=======================================================================
% spm_get('Add',H)
%fprintf('Doing: spm_get(''Add''...\n')
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
	H=P2;
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
	FPath  = [setstr(ones(nItems,1)*[WDir,'/']),Items];
	if isempty(P), P=FPath; else, P=str2mat(P,FPath); end
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
% (The str2mat bit ensures single strings are handled by column)
if size(P,1), P(:,all(str2mat(P,' ')==' '))=[]; end

%-Return new P to holding object
%-----------------------------------------------------------------------
set(findobj(F,'Tag','P'),'UserData',P);

%-Update StatusLine
%-----------------------------------------------------------------------
spm_get('StatusLine',size(P,1),n)

return


elseif strcmp(Action,'Delete')
%=======================================================================
% spm_get('Delete',h)
%fprintf('Doing: spm_get(''Delete'')\n')
F      = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
if nargin<2, h=gcbo; else, h=P2; end
WDir   = get(findobj(F,'Tag','WDir'),'UserData');
P      = get(findobj(F,'Tag','P'),'UserData');
n      = get(findobj(F,'Tag','Prompt'),'UserData');
IName  = get(h,'String');
IName(1:find(IName==':')) = [];
Items  = get(h,'UserData');
nItems = size(Items,1);
FPath  = [setstr(ones(nItems,1)*[WDir,'/']),Items];
nP     = size(P,1);

%-Compare Items with end of P (Can only delete from the end)
tmp    = str2mat(P(nP-(nItems-1):nP,:),FPath);
if all(all(tmp(1:nItems,:)==tmp(nItems+1:2*nItems,:)))
	P(nP-(nItems-1):nP,:)=[];
	if ~isempty(P)
		%-Chop off any redundant trailing spaces
		P(:,all(str2mat(P,' ')==' '))=[];
	end
	set(findobj(F,'Tag','P'),'UserData',P);
	set(h,'String',IName,'Color','k',...
		'Tag','IName',...
		'ButtonDownFcn','spm_get(''Add'')')
	spm_get('StatusLine',size(P,1),n)
end
return


elseif strcmp(Action,'Reset')
%=======================================================================
% spm_get('Reset')
%fprintf('Doing: spm_get(''Reset'')\n')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(findobj(F,'Tag','P'),'UserData',[]);
spm_get('dir')
spm_get('StatusLine')
return


elseif strcmp(Action,'All')
%=======================================================================
% spm_get('All')
%fprintf('Doing: spm_get(''All'')\n')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(F,'Pointer','Watch')
H = flipud(findobj(get(gca,'Children'),'Flat','Tag','IName'));
spm_get('Add',H);
set(F,'Pointer','Arrow')
return


elseif strcmp(Action,'CmdLine')
%=======================================================================
% P = spm_get('CmdLine',n,Prompt,P,WDir,AllowEnd)
%fprintf('Doing: spm_get(''CmdLine''...\n')
if nargin<6, AllowEnd=0; else AllowEnd=P6; end
if nargin<5, WDir=''; else WDir=P5; end
if isempty(WDir), WDir=pwd; else, WDir=deblank(WDir); end
if nargin<4, P=[]; else P=P4; end
if nargin<3, Prompt='Select files...'; else Prompt=P3; end
if nargin<2, n=Inf; else n=P2; end

if n==0, return, end

nP=size(P,1);

fprintf('\n'), fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
fprintf('\t%s\n',Prompt)
fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
fprintf('Current working directory: %s\n',WDir)
fprintf('(Prepended to relative paths)\n')
if nP>0
	fprintf('\nAlready selected:\n')
	for pP = [1:size(P,1)], fprintf('\t%s\n',P(pP,:)), end
end
fprintf('\nEnter paths :')
AllowEnd = AllowEnd | isinf(n);
Tstr = 'END';
if AllowEnd, fprintf(' (Type "END" to terminate input.)')
	else fprintf(' to %d items:',abs(n)-nP), end

Done=0;
while ~Done
	%-Get input
	str = []; while isempty(str)
		str=deblank(input(sprintf('  %3d  : ',nP+1),'s'));end
	%-Prepend WDir to incomplete pathnames
	if (str(1)~='/')&(~strcmp(str,Tstr)), str = [WDir,'/',str]; end
	if (n>0),OK=exist(str)==2;else,OK=~unix(['test -d ',str]);end
	while ~( OK | (strcmp(str,Tstr) & AllowEnd) )
		if strcmp(str,Tstr)
			fprintf('%c\tSelect %d files!',7,abs(n))
			else, fprintf('%c\t%s doesn''t exist!',7,str), end
		str=[]; while isempty(str)
			str=deblank(input(sprintf('  %3d  : ',nP+1),'s'));end
		if (str(1)~='/')&(~strcmp(str,Tstr)), str=[WDir,'/',str]; end
		if (n>0),OK=exist(str)==2;else,OK=~unix(['test -d ',str]);end
	end
	if ~strcmp(str,Tstr)
		if isempty(P), P=str; else, P=str2mat(P,str); end, end
	nP = nP+1;
	Done = (strcmp(str,Tstr) | (nP==abs(n)));
end % (while)

R1 = P;
return


elseif strcmp(Action,'GUI2CmdLine')
%=======================================================================
% P = spm_get('GUI2CmdLine')
%fprintf('Doing: spm_get(''GUI2CmdLine'')\n')

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
	spm_get('StatusLine',size(PNew,1),n)
end

delete(Handles);
return


elseif strcmp(Action,'Edit')
%=======================================================================
% spm_get('Edit')
%fprintf('Doing: spm(''Edit'')\n')

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
return


elseif strcmp(Action,'EditDone')
%=======================================================================
% spm_get('EditDone',OK)
%fprintf('Doing: spm_get(''EditDone''...\n')
if nargin<2, OK=0; else, OK=strcmp(P2,'OK'); end
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

EditHandles = get(findobj(F,'Tag','EditHandles'),'UserData');
set(EditHandles,'Visible','off')
h_EditWindow = findobj(F,'Tag','EditWindow');

if OK & get(h_EditWindow,'UserData')
	%-Only update if OK was pressed, and there was an edit.
	P = get(h_EditWindow,'String');
	n = get(findobj(F,'Tag','Prompt'),'UserData');
	%-Condition P, by removing blank lines and truncating to n items.
	P = P(~all(P'==0),:);
	if (finite(n))&(~isempty(P)), P=P(1:min(size(P,1),abs(n)),:); end
	%-Write P to 'P' 'Tag'ged object
	set(findobj(F,'Tag','P'),'UserData',P);
	%-Redo directory listing
	spm_get('dir')
	%-Redo status line
	spm_get('StatusLine',size(P,1),n)
end
delete(EditHandles), drawnow

return


elseif strcmp(Action,'Done')
%=======================================================================
% spm_get('Done')
%fprintf('Doing: spm_get(''Done'')\n')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
n = get(findobj(F,'Tag','Prompt'),'UserData');

%-If #files was specified, and not enough have been selected, then return
if finite(n)
	P = get(findobj(F,'Tag','P'),'UserData');
	if (size(P,1)<abs(n)), fprintf('%c',7), return, end
end

%-Done, set Done UserData tag for handling
set(findobj(F,'Tag','Done'),'UserData',1)
return


else
%=======================================================================
error('Illegal Action string')


%=======================================================================
end

