function R1 = spm_get_cb(Action,P2,P3,P4,P5,P6)
% Callback and setup function for spm_get.
% FORMAT R1 = spm_get_cb(Action,P2,P3,P4,P5,P6)
%           - A multi function function!
%
% FORMAT F = spm_get_cb('CreateFig',LastDirs,Filter)
% Sets up the file selection window in next avaliable figure.
% F        - Figure used.
% LastDirs - String matrix for directory history.
% Filter   - Filename filter.
%
% FORMAT F = spm_get_cb('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.
% F        - Figure used.
% Vis      - Figure visibility (string, 'on' or 'off')
% n        - Number of filepaths/directories to obtain
% Prompt   - Prompt string.
% Filter   - Filename filter.
% WDir     - [optional] New Working directory
%
% FORMAT spm_get_cb('StatusLine',nP,n)
% Updates FileSelWin status line
% nP       - #Items already selected
% n        - #Items to select
%
% FORMAT spm_get_cb('cd',NewDir,bNoDir)
% Callback function for change directory objects
% NewDir   - New Working directory (if not held in current object)
% bNoDir   - Supresses listing of new directory if true
%
% FORMAT spm_get_cb('dir',WDir)
% Lists current working directory and attatches callbacks to contents
% WDir     - Working directory (Defaults to UserData of WDir Tagged object.
%
% FORMAT spm_get_cb('Add',h)
% Adds current object (of object h) to selection
% h        - Object handle of an item
%
% FORMAT spm_get_cb('Delete')
% Deletes current object from the selection, if it was the last chosen
%
% FORMAT spm_get_cb('Reset')
% Resets the selection and redisplays the directory
%
% FORMAT spm_get_cb('All')
% Adds all remaining unselected items to the selection
%
% FORMAT P = spm_get_cb('CmdLine',n,Prompt,P,WDir)
% Uses command line to get items
% n        - # items to get
% Prompt   - Prompt string
% P        - Currently selected items
% WDir     - Working directory to use (defaults to pwd)
%
% FORMAT P = spm_get_cb('GUI2CmdLine')
% GUI gateway to command line input.
% P        - Currently selected items
%
% FORMAT spm_get_cb('Edit')
% Invokes window for editing currently selected items
%
% FORMAT spm_get_cb('EditDone',OK)
% Processes the results of the edit.
% OK       - 'OK' or 'CANCEL' - Edit status.
%
% FORMAT spm_get_cb('Done')
% Callback for "Done" button.
%
%_______________________________________________________________________
%
% spm_get_cb is the "engine" for spm_get.
%
% This function requires spm_list_files, and replaces SPM94 functions
% spm_all and spm_ls, which are now obsolete.
%
% See also spm_get
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%



if strcmp(Action,'CreateFig')
%=======================================================================
% F = spm_get_cb('CreateFig',LastDirs,Filter)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<3 Filter=[]; else Filter=P3; end
if nargin<2 LastDirs=[]; else LastDirs=P2; end
if isempty(LastDirs), LastDirs=pwd; end
LastDirs=str2mat(LastDirs,getenv('HOME'));
if (exist('spm_get_path')==2), LastDirs=str2mat(LastDirs,spm_get_path); end
	

%-Create window, compute scaling for screen size
%-----------------------------------------------------------------------
S = get(0,'ScreenSize');
A = [S(3)/1152 S(4)/900 S(3)/1152 S(4)/900];
F = figure('Name','','Tag','SelFileWin',...
	'NumberTitle','off',...
	'Color',[1 1 1]*.8,...
	'Resize','off',...
	'Visible','off',...
	'Units','Pixels',...
	'Position',[108 008 400 395].*A);
R1=F;

%-User control objects with callbacks
%-----------------------------------------------------------------------

uicontrol(F,'Style','Frame','Tag','P','UserData',[],...
	'Position',[001 271 400 125].*A);

uicontrol(F,'Style','Frame',...
	'BackgroundColor',[1 1 1]*.8,...
	'Position',[010 355 380 030].*A);

uicontrol(F,'Style','Text','Tag','Prompt','UserData',[],...
	'String','<Prompt not set yet>',...
	'ForegroundColor','k',...
	'BackgroundColor',[1 1 1]*.8,...
	'HorizontalAlignment','Center',...
	'Position',[020 360 360 020].*A);

WDir=deblank(LastDirs(1,:));
uicontrol(F,'Style','Edit','String',WDir,...
	'Tag','WDir','UserData',WDir,...
	'ForegroundColor','r','BackgroundColor',[.8,.8,1],...
	'HorizontalAlignment','Left',...
	'CallBack','spm_get_cb(''cd'')',...
	'Position',[010 330 380 020].*A);

uicontrol(F,'Style','PopUp','Tag','LastDirsPopup',...
	'String',str2mat('Previous Directories...',LastDirs),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','r',...
	'UserData',size(LastDirs,1),...
	'Callback','spm_get_cb(''cd'')',...
	'Position',[010-1 310-1 380+1 020].*A);

%	'Callback','spm_get_cb(''cd'',get(gco,''Value''))',...

%uicontrol(F,'Style','PopUp','Tag','SubDirsPopup',...
%	'HorizontalAlignment','Center',...
%	'ForegroundColor','r',...
%	'String','subdirectories...',...
%	'Callback','spm_get_cb(''cd'',get(gco,''Value''))',...
%	'Position',[020-1 290 360+1 020].*A);

uicontrol(F,'Style','text','String','Filter:',...
	'Position',[010 280 040 020].*A);

uicontrol(F,'Style','Edit','Tag','Filter',...
	'String','*',...
	'UserData','*',...
	'BackgroundColor',[.8,.8,1],...
	'Callback','spm_get_cb(''dir'')',...
	'Position',[050 280 095 020].*A);

uicontrol(F,'Style','Pushbutton','String','All',...
	'Callback','spm_get_cb(''All'')',...
	'Position',[150 280 030 020].*A);

uicontrol(F,'Style','Pushbutton','String','Edit',...
	'ForegroundColor','b',...
	'Callback','spm_get_cb(''Edit'')',...
	'Position',[185 280 040 020].*A);

uicontrol(F,'Style','Pushbutton','String','Keybd',...
	'ForegroundColor','k',...
	'Callback','spm_get_cb(''GUI2CmdLine'')',...
	'Position',[230 280 050 020].*A);

uicontrol(F,'Style','Pushbutton','String','Reset',...
	'ForegroundColor','r',...
	'Callback','spm_get_cb(''Reset'')',...
	'Position',[285 280 050 020].*A);

uicontrol(F,'Style','Pushbutton','String','Done',...
	'ForegroundColor','m',...
	'Tag','Done','UserData',0,...
	'Callback','spm_get_cb(''Done'')',...
	'Position',[340 280 050 020].*A);

uicontrol(F,'Style','Pushbutton','String','/\',...
	'Callback','spm_position(gca,1,[0 -0.5]);',...
	'Position',[370 245 020 020].*A);

uicontrol(F,'Style','Pushbutton','String','\/',...
	'Callback','spm_position(gca,1,[0 0.5]);',...
	'Position',[370 040 020 020].*A);

uicontrol(F,'Style','Frame','Tag','StatusArea',...
	'Position',[001 001 400 030].*A);

uicontrol(F,'Style','Text','Tag','StatusLine',...
	'String','<Not set yet>',...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[020 005 360 020].*A);

return

elseif strcmp(Action,'Initialise')
%=======================================================================
% F = spm_get_cb('Initialise',Vis,n,Prompt,Filter,WDir)
% (Re)Initialise SelFileWin, create one if necessary.

if nargin<6 WDir=''; else WDir=P6; end
if nargin<5 Filter=0; else Filter=P5; end
if nargin<4 Prompt='Select files...'; else Prompt=P4; end
if nargin<3 n=Inf; else n=P3; end
if nargin<2 Vis='on'; else Vis=P2; end

%-Recover SelFileWin figure number
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');

%-If closing SelFileWin then close and return
if strcmp('Vis','close'), if ~isempty(F), close(F), end, R1=[]; return, end

%-If no SelFileWin, then create one
if isempty(F), F=spm_get_cb('CreateFig'); end
figure(F)

%-If Vis=='off', make invisible and delete 'dir' text axes
if strcmp(Vis,'off'), set(F,'Visible','off'), delete(gca), end

%-clear P
set(findobj(F,'Tag','P'),'UserData',[])

%-Reset prompt, & n, stored in UserData
set(findobj(F,'Tag','Prompt'),'String',Prompt,'UserData',n)

%-Reset StatusLine
spm_get_cb('StatusLine',0,n)

%-Done UserData to zero
set(findobj(F,'Tag','Done'),'UserData',0)

%-Filter string
%-Only set if different to previously passed one, held as UserData of
% Filter Edit uicontrol.
if isstr(Filter)
	if ~strcmp(Filter,get(findobj(F,'Tag','Filter'),'UserData'))
		set(findobj(F,'Tag','Filter'),'String',Filter,'UserData',Filter)
	end
end

%-Change to new working directory if one was specified
if ~isempty(WDir), spm_get_cb('cd',WDir,1), end

%-Make visible and do 'dir', if Vis='on'
if strcmp(Vis,'on'), spm_get_cb('dir'), set(F,'Visible','on'), end

%-Return figure handle
R1 = F;

return

elseif strcmp(Action,'StatusLine')
%=======================================================================
% spm_get_cb('StatusLine',nP,n)
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
%-P2 = Specific directory
%-P3 = Boolean, supresses directory listing if true.
%-Condition arguments
%-----------------------------------------------------------------------
if (nargin<3), bDirList=1; else, bDirList = ~P3; end
if (nargin<2), NewDir=''; else, NewDir=deblank(P2); end

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(F,'Pointer','Watch')

if isempty(NewDir)
	%-Current object contains new directory, either UIEdit, UIPopMenu or text
	NewDir=get(gco,'String');
	if strcmp(get(gco,'Tag'),'LastDirsPopup')
		val = get(gco,'Value');
		if val==1 set(F,'Pointer','Arrow'), return, end
		NewDir=NewDir(val,:);
	end
	NewDir=deblank(NewDir);
end % (if)

%-Condition relative pathnames
%-----------------------------------------------------------------------
if NewDir(1)~='/'
	WDir=get(findobj(F,'Tag','WDir'),'UserData');
	if NewDir(1)~='.'
		%-Prepend current working directory
		if strcmp(WDir,'/') NewDir=['/',NewDir];
			else NewDir=[WDir,'/',NewDir]; end
	else
		%-Sort out "dot" options...
		NewDir(1)=[];
		if isempty(NewDir), set(F,'Pointer','Arrow'), return, end
		if NewDir(1)=='.'
			%-NewDir started "..", go down a directory from WDir
			WDir=WDir(1:max(find(WDir=='/'))-1);
			if isempty(WDir), WDir='/'; end
			NewDir(1)=[];
		end
		NewDir=[WDir,NewDir];
	end % (if NewDir(1)~='.')
end % (if NewDir...)

%-Changing directory, delete current directory listing (Done again in
% Action=='dir', but neater to do here in advance on this occasion)
figure(F), delete(gca)


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
	LastDirs(IDRows,:)=[];
end % (if)

%-Delate NewDir from top of LastDirs if present in fixed part of LastDirs
if any( all(FLastDirs'==setstr(ones(size(FLastDirs,1),1)*LastDirs(1,:))') )
	LastDirs(1,:)=[]; end

%-Lose directories from the bottom of variable LastDirs if oversize
Vs = 10; %-Maximum allowable number of variable LastDirs entries
if size(LastDirs,1)>Vs LastDirs(Vs+1:size(LastDirs,1),:)=[]; end

%-Put LastDirs and FLastDirs back together
LastDirs = [LastDirs; FLastDirs];


%-Write LastDirs to LastDirsPopup, write new directory to WDir Edit widget
%-----------------------------------------------------------------------
set(findobj(F,'Tag','LastDirsPopup'),...
	'String',str2mat('Previous Directories...',LastDirs),...
	'Value',1)

set(findobj(F,'Tag','WDir'),'String',NewDir,'UserData',NewDir)


%-Display new directory, unless supressed
%-----------------------------------------------------------------------
if bDirList, spm_get_cb('dir'), end
	
drawnow, set(F,'Pointer','Arrow')

return

elseif strcmp(Action,'dir')
%=======================================================================
% spm_get_cb('dir',WDir)
%
% Creates list of text objects with an associated
% ButtonDownFcn functions that call spm_get_cb for processing
%
%-WDir - Directory, defaults to UserData of 'WDir' tagged object in
% SelFileWin figure

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
if nargin<2, WDir=get(findobj(F,'Tag','WDir'),'UserData');
	else, WDir = P2; end

n      = get(findobj(F,'Tag','Prompt'),'UserData');
Filter = get(findobj(F,'Tag','Filter'),'String');
if isempty(Filter)
	Filter='*'; set(findobj(F,'Tag','Filter'),'String',Filter), end

%-Write current directory to WDir widget, if new directory specified
%-----------------------------------------------------------------------
if nargin>1
	set(findobj(F,'Tag','WDir'),'String',WDir,'UserData',WDir), end

%-Set up an axis and get the list of filenames
%-----------------------------------------------------------------------
figure(F)
delete(gca)
axes('Position',[0 0 1 1],'Visible','off')
[Files,Dirs] = spm_list_files(WDir,Filter);
q      = .64;					% y position of pathnames
dy     = 0.046;					% spacing of text objects

%-List current directory
%-----------------------------------------------------------------------
if isempty(Dirs)
	text(0.01,q,'Permission denied, or non-existent directory',...
		'FontWeight','bold','Color','r');
end

%-Create list of directories with appropriate 'ButtonDownFcn': -
%-----------------------------------------------------------------------
y     = q;
for i = 1:size(Dirs,1)
	text(0.01,y,Dirs(i,:),'Tag','DirName',...
		'FontWeight','bold','Color','r',...
		'ButtonDownFcn','spm_get_cb(''cd'')');
	y = y - dy;
end

%-Create list of directories in pulldown menu
%-----------------------------------------------------------------------
%set(findobj(F,'Tag','SubDirsPopup'),'String',...
%	str2mat('SubDirectories...',Dirs))


Items=Files; if (n<0) Items=Dirs; end

%-Create list of filenames with appropriate 'ButtonDownFcn': -
%-----------------------------------------------------------------------
y     = q;
for i = 1:size(Items,1)
	text(0.40,y,Items(i,:),'Tag','IName',...
		'Color','k',...
		'ButtonDownFcn','spm_get_cb(''Add'')');
	y = y - dy;
end


return


elseif strcmp(Action,'Add')
%=======================================================================
% spm_get_cb('Add',h)
% Add filename to current list
%-h - [Optional] handle to object
%
if (nargin==2), h=P2; else h=gco; end

%-Recover variables from holding objects
%-----------------------------------------------------------------------
F    = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
n    = get(findobj(F,'Tag','Prompt'),'UserData');
P    = get(findobj(F,'Tag','P'),'UserData');
WDir = get(findobj(F,'Tag','WDir'),'UserData');

%-If #items has been specified, and P is full, then return
%-----------------------------------------------------------------------
if finite(n), if (size(P,1)>=abs(n)), return, end, end

%-Add item to P, update StatusLine
%-----------------------------------------------------------------------
FName = get(h,'String');
FPath = [WDir,'/',FName];
if isempty(P) P=FPath; else P=str2mat(P,FPath); end

set(findobj(F,'Tag','P'),'UserData',P);

set(h,'String',[int2str(size(P,1)),' :',FName],'Color','b',...
	'Tag','SelIName',...
	'ButtonDownFcn','spm_get_cb(''Delete'')')

spm_get_cb('StatusLine',size(P,1),n)


return


elseif strcmp(Action,'Delete')
%=======================================================================
% spm_get_cb('Delete')
F     = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
WDir  = get(findobj(F,'Tag','WDir'),'UserData');
P     = get(findobj(F,'Tag','P'),'UserData');
n     = get(findobj(F,'Tag','Prompt'),'UserData');
FName = get(gco,'String'); FName(1:find(FName==':')) = [];
FPath = [WDir,'/',FName];
NoP   = size(P,1);
if strcmp(P(NoP,:),FPath)
	P(NoP,:)=[];
	set(findobj(F,'Tag','P'),'UserData',P);
	set(gco,'String',FName,'Color','k',...
		'Tag','IName',...
		'ButtonDownFcn','spm_get_cb(''Add'')')
	spm_get_cb('StatusLine',size(P,1),n)
end


elseif strcmp(Action,'Reset')
%=======================================================================
% spm_get_cb('Reset')
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
set(findobj(F,'Tag','P'),'UserData',[]);
spm_get_cb('dir')
spm_get_cb('StatusLine')

return

elseif strcmp(Action,'All')
%=======================================================================
% spm_get_cb('All')
%
%-MatLab's built in findobj gives a bus error for more than 99 matching
% target objects. :-(
h = get(gca,'Children');
Handles = [];
while length(h)>99
	Handles = [Handles; findobj(h(1:99),'Flat','Tag','IName')];
	h(1:99)=[];
end
Handles = flipud([Handles; findobj(h,'Flat','Tag','IName')]);

for h = Handles'
	spm_get_cb('Add',h);
end

return

elseif strcmp(Action,'CmdLine')
%=======================================================================
% P = spm_get_cb('CmdLine',n,Prompt,P,WDir,AllowEnd)
if nargin<6, AllowEnd=0; else AllowEnd=P6; end
if nargin<5, WDir=pwd; else WDir=P5; end
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
	str = [input(sprintf('  %3d  : ',nP+1),'s'),' '];
	%-Prepend WDir to incomplete pathnames, prevents exist searching the path.
	if (str(1)~='/')&(~strcmp(str,Tstr)), str = [WDir,'/',str]; end
	if (n>0), OK=exist(str)==2; else, OK=~unix(['test -d ',str]); end
	while ~( OK | (strcmp(str,Tstr) & AllowEnd) )
		fprintf('%c\t%s doesn''t exist!',7,str)
		str = [input(sprintf('  %3d! : ',nP+1),'s'),' '];
		if (str(1)~='/')&(~strcmp(str,Tstr)), str = [WDir,'/',str]; end
		if (n>0), OK=exist(str)==2; else, OK=~unix(['test -d ',str]); end
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
% P = spm_get_cb('GUI2CmdLine')

F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
Prompt = get(findobj(F,'Tag','Prompt'),'String');
n      = get(findobj(F,'Tag','Prompt'),'UserData');
WDir   = get(findobj(F,'Tag','WDir'),'UserData');
P      = get(findobj(F,'Tag','P'),'UserData');

A = get(F,'Position'); A = [A(3)/400,A(4)/395,A(3)/400,A(4)/395];
Handles = [];
h = uicontrol(F,'Style','Frame',...
		'Position',[000 000 400 395].*A);
Handles = [Handles, h];

h = uicontrol(F,'Style','Frame',...
	'BackgroundColor',[1 1 1]*.8,...
	'Position',[010 355 380 030].*A);
Handles = [Handles, h];

h = uicontrol(F,'Style','Text','String',Prompt,...
	'ForegroundColor','k',...
	'BackgroundColor',[1 1 1]*.8,...
	'HorizontalAlignment','Center',...
	'Position',[020 360 360 020].*A);
Handles = [Handles, h];

h = uicontrol(F,'Style','Frame',...
	'Position',[001 001 400 030].*A);
Handles = [Handles, h];

h = uicontrol(F,'Style','Text',...
	'String','Use command window...',...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[020 005 360 020].*A);
Handles = [Handles, h];

PNew = spm_get_cb('CmdLine',n,Prompt,P,WDir,1);

if ~strcmp(P(:),PNew(:))
	set(findobj(F,'Tag','P'),'UserData',PNew);
	spm_get_cb('dir')
	spm_get_cb('StatusLine',size(PNew,1),n)
end

delete(Handles);

return

elseif strcmp(Action,'Edit')
%=======================================================================
F = findobj(get(0,'Children'),'Flat','Tag','SelFileWin');
n = get(findobj(F,'Tag','Prompt'),'UserData');

set(F,'Units','Pixels');
A = get(F,'Position'); A = [A(3)/400,A(4)/395,A(3)/400,A(4)/395];

EditHandles = [];

h = uicontrol(F,'Style','Frame',...
	'Tag','EditHandles',...
	'Position',[001 001 400 395].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Frame',...
	'BackgroundColor',[1 1 1]*.8,...
	'Position',[010 355 380 030].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Text',...
	'String',get(findobj(F,'Tag','Prompt'),'String'),...
	'ForegroundColor','k',...
	'BackgroundColor',[1 1 1]*.8,...
	'HorizontalAlignment','Center',...
	'Position',[020 360 360 020].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Text',...
	'String','Edit the filename window...',...
	'ForegroundColor','w',...
	'HorizontalAlignment','Left',...
	'Position',[010 325 210 020].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Pushbutton','String','Cancel',...
	'ForegroundColor','r',...
	'Callback','spm_get_cb(''EditDone'',''Cancel'')',...
	'Position',[220 325 070 020].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Pushbutton','String','OK',...
	'ForegroundColor','m',...
	'Callback','spm_get_cb(''EditDone'',''OK'')',...
	'Position',[310 325 070 020].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Edit',...
	'String',get(findobj(F,'Tag','P'),'UserData'),...
	'Tag','EditWindow','UserData',0,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','b','BackgroundColor',[.8,.8,1],...
	'Max',2,...
	'Callback','set(gco,''UserData'',1)',...
	'Position',[010 040 380 275].*A);
EditHandles = [EditHandles, h];

h = uicontrol(F,'Style','Frame',...
	'Position',[001 001 400 030].*A);
EditHandles = [EditHandles, h];

str=[]; if finite(n), str=['Use only top ',int2str(abs(n)),' lines.  ']; end
h = uicontrol(F,'Style','Text',...
	'String',[str,'(Blank lines will be removed)'],...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w',...
	'Position',[020 005 360 020].*A);
EditHandles = [EditHandles, h];

set(findobj(F,'Tag','EditHandles'),'UserData',EditHandles)

return

elseif strcmp(Action,'EditDone')
%=======================================================================
% spm_get_cb('EditDone',OK)
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
	spm_get_cb('dir')
	%-Redo status line
	spm_get_cb('StatusLine',size(P,1),n)
end
delete(EditHandles), drawnow


return


elseif strcmp(Action,'Done')
%=======================================================================
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
