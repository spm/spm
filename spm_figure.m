function R1=spm_figure(Action,P2,P3,P4,P5)
% Setup and callback functions for Graphics window
% FORMAT R1=spm_figure(Action,P2,P3,P4,P5)
%	- An embedded callback, multi-function function
%       - For detailed programmers comments, see the format specifications
%         below
%_______________________________________________________________________
%
% spm_figure creates and manages the 'Graphics' window. This window and
% these facilities may be used independently of SPM, and any number of
% Graphics windows my be used within the same MatLab session. (Though
% only one SPM 'Graphics' 'Tag'ed window is permitted.
%
% The Graphics window is provided with a menu bar at the top that
% facilitates editing and printing of the current graphic display,
% enabling interactive editing of graphic output prior to printing
% (e.g. selection of color maps, deleting, moving and editing graphics
% objects or adding text)
%
% Print: Creates a footnote (detailing the SPM version, user & date)
% and evaluates PRINTSTR (see spm_defaults.m). Graphics windows with
% multi-page axes are printed page by page.
%
% Clear: Clears the Graphics window, leaving the menu bar. If in SPM
% usage (figure 'Tag'ed as 'Graphics' then the 'Interactive' 'Tag'ed
% window is also cleared, and it's name and pointer reset.
%
% Colormap options:
% * Gray and Hot: Sets the colormap to its default values and loads
%                 either a grayscale or 'hot metal' color map.
% * Split: Loads a 'split' color map from 'Split.mat' {128 x 3
%          matrix}.  The lower half is a gray scale and the upper half
%          is 'hot metal'.  This color map is used for viewing 'rendered' 
%          SPM{Z} on a PET, MRI or other background images
%
% Colormap effects:
% * Invert: Inverts (flips) the current color map.
% * Brighten and Darken: Brighten and Darken the current colourmap
% 	   using the MatLab BRIGHTEN command, with  beta's of +0.2 and -0.2
% 	   respectively.
%
% Editing: Right button ('alt' button) cancels operations
% * Cut  : Deletes the graphics object next selected (if deletable)
% * Move : To re-position a text, uicontrol or axis object using a
%          'drag and drop' implementation (i.e. depress - move - release)
%          Using the middle 'extend' mouse button on a text object moves
%          the axes containing the text - i.e. blocks of text.
% * Size : Re-sizes the text, uicontrol or axis object next selected
%          {left button - decrease, middle button  - increase} by a factor
%          of 1.24 (or increase/decrease FontSize by 2 dpi)
% * Text : Creates an editable text widget that produces a text object as
%          its CallBack.
%          This text object can then be manipulated using the edit facilities.
% * Edit : To edit text, select a text object with the circle cursor,
%          and edit the text in the editable text widget that appears.
%
%
% Objects with the attribute Tag = 'NoDelete' are exempt from deletion
% when spm_clf (spm_figure('Clear')) is used.
%
% For SPM usage, the figure should be 'Tag'ed as 'Graphics'.
%
% For SPM power users, and programmers, spm_figure provides utility
% routines for using the SPM graphics interface. Of particular use are
% the FindWin and Clear functions See the embedded callback reference
% in the main body of spm_figure, below the help text.
%
% See also: spm_print, spm_clf
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `spm_figure Create` is           )
%( equivalent to `spm('Create')`.                                      )
%
% FORMAT F = spm_figure
% [ShortCut] Defaults to Action 'Create'
%
% FORMAT F = spm_figure(F) - numeric F
% [ShortCut] Defaults to spm_figure('CreateBar',F)
%
% FORMAT F = spm_figure('Create',Tag,Name,Visible)
% Create a full length WhiteBg figure 'Tag'ed Tag (if specified),
% with a ToolBar.
% Equivalent to spm_figure('CreateWin','Tag') and spm_figure('CreateBar')
% Tag	  - 'Tag' string for figure.
% Name    - Name for window
% Visible - 'on' or 'off'
% F	  - Figure used
%
% FORMAT F = spm_figure('FindWin',F)
% Finds window with 'Tag' or figure numnber F - returns empty F if not found
% F	- (Input)  Figure to use [Optional] - 'Tag' string or figure number.
%	- Defaults to 'Graphics'
% F	- (Output) Figure number (if found) or empty (if not).
%
% FORMAT F = spm_figure('FindParentWin',h)
% Finds window containing the object whose handle is specified
% h	- Handle of object whose parent figure is required
%	- If a vector, then first object handle is used
% F	- Number or parent figure
%
% FORMAT spm_figure('Clear',F)
% F	- 'Tag' string or figure number of figure to clear, defaults to gcf
% Clears figure, leaving ToolBar (& other objects 'Tag'ed as 'NoDelete')
% intact. If figure F is 'Tag'ged 'Interactive' (SPM usage), then the window
% name and pointer are reset.
%
% FORMAT F = spm_figure('CreateWin',Tag,Name,Visible)
% Creates a full length WhiteBg figure 'Tag'ged Tag (if specified).
% F	  - Figure created
% Tag	  - Tag for window
% Name    - Name for window
% Visible - 'on' or 'off'
%
% FORMAT spm_figure('CreateBar',F)
% Creates toolbar in figure F (defaults to gcf). F can be a 'Tag'
% If the figure is 'Tag'ed as 'Graphics' (SPM usage), then the Print button
% callback is set to attempt to clear an 'Interactive' figure too.
%
% FORMAT spm_figure('Print',F,PFile)
% SPM print function: Appends footnote & executes PRINTSTR
% F	- [Optional] Figure to print. ('Tag' or figure number)
%	  Defaults to figure 'Tag'ed as 'Graphics'.
%	  If none found, uses CurrentFigure if avaliable.
% PFile - [Optional] File to print to.
%	  If specified then PRINTSTR is overridden.
% PRINTSTR - global variable holding print command to be evaluated
%	  Defaults to 'print -dps2 fig.ps'
% If objects 'Tag'ed 'NextPage' and 'PrevPage' are found, then the
% pages are shown and printed in order. In breif, pages are held as
% seperate axes, with ony one 'Visible' at any one time. The handles of
% the "page" axes are stored in the 'UserData' of the 'NextPage'
% object, while the 'PrevPage' object holds the current page number.
% See spm_help('!Disp') for details on setting up paging axes.
%
% FORMAT spm_figure('NewPage',hPage,MoveOn)
% SPM pagination function: Makes objects with handles hPage paginated
% Creates pagination buttons if necessary.
% hPage	 - Handles of objects to stick to this page
% MoveOn - If specified, causes hPage objects to be made invisible
%
% FORMAT spm_figure('DeletePageControls',F)
% SPM pagination function: Deletes page controls
% F	- [Optional] Figure in which to attempt to turn the page
%         Defaults to 'Graphics' 'Tag'ged window
%
% FORMAT spm_figure('WaterMark',F,str,Tag,Angle)
% Adds watermark to figure windows.
% F	- Figure for watermark. Defaults to gcf
% str   - Watermark string. Defaults (missing or empty) to SPM
% Tag   - Tag for watermark axes. Defaults to ''
% Angle - Angle for watermark. Defaults to -45
%
% FORMAT spm_figure('TurnPage',move,F)
% SPM pagination function: Turn to specified page
%
% FORMAT spm_figure('ColorMap')
% Callback for "ColorMap" buttons
%
% FORMAT spm_figure('GraphicsCut')
% Callback for "Cut" button
%
% FORMAT spm_figure('GraphicsMoveStart')
% Callback for "Move" button
%
% FORMAT spm_figure('GraphicsMoveMotion')
% Callback for move function
%
% FORMAT spm_figure('GraphicsMoveEnd')
% Callback for move function
%
% FORMAT spm_figure('GraphicsSize')
% Callback for "Size" button
%
% FORMAT spm_figure('GraphicsText')
% Callback for "Text" button
%
% FORMAT spm_figure('GraphicsTextEdit')
% Callback for "Edit" button
%_______________________________________________________________________


%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'Create'; nargin=1; end
if ~isstr(Action), P2=Action; Action='CreateBar'; end

if strcmp(lower(Action),lower('Create'))
%=======================================================================
% F = spm_figure('Create',Tag,Name,Visible)
%-Condition arguments
if nargin<4, Visible='on'; else, Visible=P4; end
if nargin<3, Name='Graphics'; else, Name=P3; end
if nargin<2, Tag = []; else, Tag=P2; end

F = spm_figure('CreateWin',Tag,Name,Visible);
spm_figure('CreateBar',F)
R1 = F;
return


elseif strcmp(lower(Action),lower('CreateWin'))
%=======================================================================
% F=spm_figure('CreateWin',Tag,Name,Visible)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<4, Visible=''; else, Visible=P4; end
if isempty(Visible), Visible='on'; end
if nargin<3, Name='Graphics'; else, Name=P3; end
if nargin<2, Tag=[]; else, Tag=P2; end

S0   = get(0,'ScreenSize');
WS   = [S0(3)/1152 S0(4)/900 S0(3)/1152 S0(4)/900];

F      = figure('Position',[515 008 600 865].*WS,...
	'Resize','off',...
	'Visible','off',...
	'PaperPosition',[.75 1.5 7 9.5],...
	'PaperType','a4letter',...
	'DefaultTextFontSize',2*round(12*min(WS)/2),...
	'Color','w',...
	'ColorMap',gray(64));
whitebg(F,'w')
if ~isempty(Tag), set(F,'Tag',Tag), end
if ~isempty(Name), set(F,'Name',Name,'NumberTitle','off'), end
set(F,'Visible',Visible)
R1 = F;
return


elseif strcmp(lower(Action),lower('FindWin'))
%=======================================================================
% F=spm_figure('FindWin',F)
% F=spm_figure('FindWin',Tag)
%-Find window: Find window with 'Tag' attribute / FigureNumber#
%-Returns empty if window cannot be found - deletes multiple tagged figs.

if nargin<2, F='Graphics'; else, F=P2; end

if isempty(F)
	% Leave F empty
elseif isstr(F)
	% Finds Graphics window with 'Tag' string - delete multiples
	Tag=F;
	F = findobj(get(0,'Children'),'Flat','Tag',Tag);
	if ( length(F) > 1 )
		%-Multiple Graphics windows - close all but most recent
		tmp   = F(2:length(F));
		close(tmp)
		F = F(1);
	end
else
	% F is supposed to be a figure number - check it
	if ~any(F==get(0,'Children')), F=[]; end
end

R1 = F;
return

elseif strcmp(lower(Action),lower('FindParentWin'))
%=======================================================================
% F=spm_figure('FindParentWin',h)
if nargin<2, error('No object specified'), else, h=P2; end
F = get(h(1),'Parent');
while ~strcmp(get(F,'Type'),'figure'), F=get(F,'Parent'); end
R1 = F;
return


elseif strcmp(lower(Action),lower('Clear'))
%=======================================================================
% spm_figure('Clear',F)
%-Clear window, leaving 'NoDelete' 'Tag'ed objects, reset pointer & name

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=P2; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Clear figure, leaving 'NoDelete' 'Tag'ed objects
%-----------------------------------------------------------------------
for h = get(F,'Children')'
	if ~strcmp(get(h,'Tag'),'NoDelete'), delete(h), end
end
set(F,'Pointer','Arrow')

%-If this is the 'Interactive' window, reset the name and pointer
if strcmp(get(F,'Tag'),'Interactive'), set(F,'Pointer','Arrow','Name',''), end

return


elseif strcmp(lower(Action),lower('Print'))
%=======================================================================
% spm_figure('Print',F,PFile)

%-Arguments & defaults
if nargin<3, PFile='spm.ps'; else, PFile=P3; end
if nargin<2, F='Graphics'; else, F=P2; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
F=spm_figure('FindWin',F);
if isempty(F), if any(get(0,'Children')), F = gcf; end, end
if isempty(F), return, end

%-See if window has paging controls
hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
iPaged    = ~isempty(hNextPage);
figure(F)

%-Retrieve print command
%-----------------------------------------------------------------------
global PRINTSTR
if isempty(PRINTSTR) | (nargin>3)
	if iPaged
		PrintCmd = ['print -dpsc2 -append ',PFile];
	else
		PrintCmd = ['print -dpsc2 ',PFile];
	end
else
	PrintCmd = PRINTSTR;
end
%-Handle specification of figure not supported my MatLab yet.
% PrintCmd = strrep(PrintCmd,'print','print -fhandle F')

%-Create footnote with SPM version, username, date and time.
%-----------------------------------------------------------------------
User  = getenv('USER');
tmp   = clock;
if exist('spm.m')==2, SPMver=spm('ver'); else, SPMver='SPM'; end
FNote = sprintf('%s (%s) - %02d/%02d/%4d (%02d:%02d)',...
		SPMver,User,tmp(3),tmp(2),tmp(1),tmp(4),tmp(5) );
%-Delete old tag lines, and print new one
delete(findobj(F,'Tag','SPMprintFootnote'));
axes('Position',[0.005,0.005,0.1,0.1],...
	'Visible','off',...
	'Tag','SPMprintFootnote')
text(0,0,FNote,'FontSize',6);

%-Print
%-----------------------------------------------------------------------
if ~iPaged
	eval(PrintCmd)
else
	hAxes     = get(hNextPage,'UserData');
	Cpage     = get(hPageNo,  'UserData');
	nPages    = size(hAxes,1);
	
	if Cpage~=1
		set(hAxes(Cpage,hAxes(Cpage,:)~=0),'Visible','off'), end
	for p = 1:nPages
		set(hAxes(p,hAxes(p,:)~=0),'Visible','on')
		eval(PrintCmd)
		set(hAxes(p,hAxes(p,:)~=0),'Visible','off')
	end
	if Cpage~=1
		set(hAxes(Cpage,hAxes(Cpage,:)~=0),'Visible','on'), end
end
return

elseif strcmp(lower(Action),lower('NewPage'))
%=======================================================================
% h = spm_figure('NewPage',hPage,MoveOn)
if nargin<3, MoveOn=0; else, MoveOn=1; end
if nargin<2, hPage=[]; else, hPage=P2(:)'; end
if isempty(hPage), error('No handles to paginate'), end

%-Work out which figure we're in
Fgraph = spm_figure('FindParentWin',hPage(1));

hNextPage = findobj(Fgraph,'Tag','NextPage');
hPrevPage = findobj(Fgraph,'Tag','PrevPage');
hPageNo   = findobj(Fgraph,'Tag','PageNo');

%-Create pagination widgets if required
%-----------------------------------------------------------------------
if isempty(hNextPage)
	hNextPage = uicontrol(Fgraph,'Style','Pushbutton',...
		'String','>',...
		'Callback','spm_figure(''TurnPage'',''+1'',gcf)',...
		'Position',[580 020 015 015].*spm('GetWinScale'),...
		'ForegroundColor',[0 0 0],...
		'Tag','NextPage','UserData',[]);
	hPrevPage = uicontrol(Fgraph,'Style','Pushbutton',...
		'String','<',...
		'Callback','spm_figure(''TurnPage'',''-1'',gcf)',...
		'Position',[565 020 015 015].*spm('GetWinScale'),...
		'Visible','on',...
		'ForegroundColor',[1 1 1]*0.5,...
		'Tag','PrevPage');
	hPageNo = uicontrol(Fgraph,'Style','Text',...
		'String','',...
		'HorizontalAlignment','center',...
		'BackgroundColor','w',...
		'Position',[560 005 040 015].*spm('GetWinScale'),...
		'Visible','on',...
		'Tag','PageNo','UserData',1);
end

%-Add handles for this page to UserData of hNextPage
%-----------------------------------------------------------------------
hAxes     = get(hNextPage,'UserData');
tmp       = max(size(hAxes,2),length(hPage));
hAxes     = [ [hAxes, zeros(size(hAxes,1),tmp-size(hAxes,2))];...
		[hPage, zeros(1,tmp-length(hPage))] ];
set(hNextPage,'UserData',hAxes)

%-Make handles for this page invisible if requested
%-----------------------------------------------------------------------
if MoveOn
	set(get(hPage,'Children'),'Visible','off')
end

%-Return handles to pagination controls if requested
if nargout>0
	R1 = [hNextPage, hPrevPage, hPageNo];
end

return


elseif strcmp(lower(Action),lower('TurnPage'))
%=======================================================================
% spm_figure('TurnPage',move,F)
if nargin<3, F='Graphics'; else, F=P3; end
if nargin<2, move=1; else, move=P2; end
Fgraph = spm_figure('FindWin',F);
if isempty(Fgraph), error('No Graphics window'), end

hNextPage = findobj(Fgraph,'Tag','NextPage');
hPrevPage = findobj(Fgraph,'Tag','PrevPage');
hPageNo   = findobj(Fgraph,'Tag','PageNo');
if isempty(hNextPage), return, end
hAxes     = get(hNextPage,'UserData');
Cpage     = get(hPageNo,  'UserData');
nPages    = size(hAxes,1);

%-Sort out new page number
if isstr(move)
	Npage = Cpage+eval(move);
else
	Npage = move;
end
Npage = max(min(Npage,nPages),1);

%-Make current page invisible, new page visible, set page number string
set(hAxes(Cpage,hAxes(Cpage,:)~=0),'Visible','off')
set(hAxes(Npage,hAxes(Npage,:)~=0),'Visible','on')
set(hPageNo,'UserData',Npage,'String',sprintf('%d / %d',Npage,nPages))

%-Set dimness of page turning controls (just for neatness)
if Npage==1, set(hPrevPage,'ForegroundColor',[1 1 1]*0.5)
else, set(hPrevPage,'ForegroundColor',[0 0 0]), end

if Npage==nPages, set(hNextPage,'ForegroundColor',[1 1 1]*0.5)
else, set(hNextPage,'ForegroundColor',[0 0 0]), end

return

elseif strcmp(lower(Action),lower('DeletePageControls'))
%=======================================================================
% spm_figure('DeletePageControls',F)
if nargin<2, F='Graphics'; else, F=P2; end
Fgraph = spm_figure('FindWin',F);
if isempty(Fgraph), error('No Graphics window'), end

hNextPage = findobj(Fgraph,'Tag','NextPage');
hPrevPage = findobj(Fgraph,'Tag','PrevPage');
hPageNo   = findobj(Fgraph,'Tag','PageNo');

delete([hNextPage hPrevPage hPageNo])

return


elseif strcmp(lower(Action),lower('WaterMark'))
%=======================================================================
% spm_figure('WaterMark',F,str,Tag,Angle)
if nargin<5, Angle=-45; else, Angle=P5; end
if nargin<4, Tag=''; else, Tag=P4; end
if nargin<3, str=''; else, str=P3; end
if isempty(str), str='SPM'; end
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=P2; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Specify watermark color from background colour
%-----------------------------------------------------------------------
Colour = get(F,'Color');
%-Only mess with grayscale backgrounds
if ~all(Colour), return, end
%-Work out colour - lighter unless grey value > 0.9
Colour = Colour+(2*(Colour(1)<0.9)-1)*0.02;

delete(findobj(F,'UserData','WaterMark'))
Units=get(F,'Units');
set(F,'Units','normalized');
figure(F)
axes('Position',[0.45,0.5,0.1,0.1],...
	'Units','normalized',...
	'Visible','off',...
	'Tag',Tag,...
	'UserData','WaterMark');
set(F,'Units',Units)
text(0.5,0.5,str,...
	'FontSize',96,...
	'FontWeight','Bold',...
	'FontName','Times',...
	'Rotation',Angle,...
	'HorizontalAlignment','Center',...
	'VerticalAlignment','middle',...
	'EraseMode','normal',...
	'Color',Colour,...
	'ButtonDownFcn',[...
		'if strcmp(get(gcf,''SelectionType''),''alt''),',...
			'delete(get(gco,''Parent'')),',...
		'end'])
return


elseif strcmp(lower(Action),lower('CreateBar'))
%=======================================================================
% spm_figure('CreateBar',F)
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=P2; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Get position and size parameters
%-----------------------------------------------------------------------
figure(F);
set(F,'Units','Pixels');
P     = get(F,'Position'); P  = P(3:4);		% Figure dimensions {pixels}
S_Gra = P./[600, 865];				% x & y scaling coefs

nBut  = 12;
nGap  = 2;
sx    = floor(P(1)./(nBut+(nGap+2)/6));		% uicontrol object width
dx    = floor(2*sx/6);				% inter-uicontrol gap
sy    = floor(20*S_Gra(1));			% uicontrol object height
x0    = dx;					% initial x position
x     = dx/2;					% uicontrol x position
y     = P(2) - sy;				% uicontrol y position
y2    = P(2) - 2.25*sy;				% uicontrol y position

%-Delete any existing 'NoDelete' ToolBar objects
%-----------------------------------------------------------------------
h = findobj(F,'Tag','NoDelete');
delete(h);

%-Create Frame for controls
%-----------------------------------------------------------------------
uicontrol(F,'Style', 'Frame',...
	'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
	'Tag','NoDelete');

%-Create uicontrol objects
%-----------------------------------------------------------------------
uicontrol(F,'String','Print' ,'Position',[x y sx sy],...
	'CallBack','spm_figure(''Print'',gcf)',...
	'Interruptible','No',...
	'Tag','NoDelete','ForegroundColor','b'); x = x+sx;

hPrint = uicontrol(F,'String','Clear' ,'Position',[x y sx sy],...
	'CallBack','spm_figure(''Clear'',gcf)',...
	'Interruptible','No',...
        'Tag','NoDelete','ForegroundColor','b'); x = x+sx+dx;
if strcmp(get(F,'Tag'),'Graphics')
	set(hPrint,'CallBack',...
	'spm_figure(''Clear'',gcf), spm_figure(''Clear'',''Interactive'')')
end

uicontrol(F,'Style','PopUp','String','ColorMap|gray|hot|split',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'spm_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'Tag','NoDelete',...
	'UserData','Pop'); x = x + 2*sx;

uicontrol(F,'Style','PopUp','String','Effects|invert|brighten|darken',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'spm_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'Tag','NoDelete',...
	'UserData','Pop'); x = x + 2*sx + dx;

uicontrol(F,'String','cut',   'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsCut'')',...
	'Interruptible','No',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','move',  'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsMoveStart'')',...
	'Interruptible','No',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','size',  'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsSize'')',...
	'Interruptible','No',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','text',  'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsText'')',...
	'Interruptible','No',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','edit',  'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsTextEdit'')',...
	'Interruptible','No',...
	'Tag','NoDelete'); x = x+sx+dx;

uicontrol(F,'String','?',     'Position',[x y dx sy],...
	'CallBack','spm_help(''spm_figure.m'')',...
	'Interruptible','No',...
	'Tag','NoDelete','ForegroundColor','g'); x = x+2*dx;
return


elseif strcmp(lower(Action),lower('ColorMap'))
%=======================================================================
% spm_figure('Colormap',ColAction,h)
if nargin<3, h=[]; else, h=P3; end
if nargin<2, ColAction='gray'; else, ColAction=P2; end
if ~isstr(ColAction)
	if ColAction==1, return, end
	Actions   = get(gco,'String');
	ColAction = deblank(Actions(ColAction,:));
end
if strcmp(ColAction,'gray')
	colormap(gray(64))
elseif strcmp(ColAction,'hot')
	colormap(hot(64))
elseif strcmp(ColAction,'split')
	load Split; colormap(split)
elseif strcmp(ColAction,'invert')
	colormap(flipud(colormap))
elseif strcmp(ColAction,'brighten')
	colormap(brighten(colormap, 0.2))
elseif strcmp(ColAction,'darken')
	colormap(brighten(colormap, -0.2))
else
	error('Illegal ColAction specification')
end

return

elseif strcmp(lower(Action),lower('GraphicsCut'))
%=======================================================================
% Delete next object clicked, provided it's deletable and in this figure
% spm_figure('GraphicsCut')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Cut: Select item to delete, RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
set(F,'Pointer','Arrow')
set(F,'Name',tmp)

if strcmp(get(gcf,'SelectionType'),'alt')
	%-Quit option
	set(hBut,'ForegroundColor','k')
	return
end

h = gco(F);
NoDel=[F; findobj(F,'Tag','NoDelete')];
if ( ~any(h==NoDel) & (gcf==F) )
	tmp  = get(h,'Type');
	if ( strcmp(tmp,'patch') | strcmp(tmp,'line') | strcmp(tmp,'image') )
		delete(get(h,'Parent'))
	else
		delete(h)
	end
end

set(hBut,'ForegroundColor','k')

return


elseif strcmp(lower(Action),lower('GraphicsMoveStart'))
%=======================================================================
% Move the next object clicked, provided it's movable and in this figure
% "alt" mouse button cancels operation
% spm_figure('GraphicsMove')

F = gcf;
hMoveBut  = gco;
set(hMoveBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Move: MiddleMouse moves blocks of text. RightMouse=cancel');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
set(F,'Name',tmp)

if strcmp(get(gcf,'SelectionType'),'alt')
	%-Quit option
	set(hMoveBut,'ForegroundColor','k')
	set(F,'Pointer','Arrow')
	return
end

hPress = gco;
NoDel=[F; findobj(F,'Tag','NoDelete')];
if (~any(hPress==NoDel) & gcf==F)
	set(F,'Units','Pixels')
	OPt  = get(F,'CurrentPoint');
	tmp  = get(hPress,'Type');
	if ( strcmp(tmp,'patch') | strcmp(tmp,'line') | strcmp(tmp,'image') )
		hMove = get(hPress,'Parent');
	elseif strcmp(get(gcf,'SelectionType'),'extend') & strcmp(tmp,'text')
		hMove = get(hPress,'Parent');
	else
		hMove = hPress;
	end
	set(hMove,'Units','Pixels');
	OPos = get(hMove,'Position');
	set(hMoveBut,'UserData',[hMove,OPt,OPos])
	set(F,'WindowButtonUpFcn','spm_figure(''GraphicsMoveEnd'')',...
		'Interruptible','no')
	set(F,'WindowButtonMotionFcn','spm_figure(''GraphicsMoveMotion'')',...
		'Interruptible','no')
	set(F,'Pointer','Fleur')
else
	set(F,'Pointer','Arrow')

end
return


elseif strcmp(lower(Action),lower('GraphicsMoveMotion'))
%=======================================================================
F     = gcf;
hPress = gco;
hMoveBut = findobj(gcf,'CallBack','spm_figure(''GraphicsMoveStart'')');
tmp   = get(hMoveBut,'UserData');
hMove = tmp(1); OPt = tmp(2:3); OPos = tmp(4:length(tmp));
CPt   = get(F,'CurrentPoint');
Disp  = CPt - OPt;
set(hMove,'Units','Pixels')
CPos  = OPos; CPos(1:2) = CPos(1:2) + Disp;
set(hMove,'Position',CPos)
set(hMove,'Units','Normalized')


elseif strcmp(lower(Action),lower('GraphicsMoveEnd'))
%=======================================================================

hMoveBut = findobj(gcf,'CallBack','spm_figure(''GraphicsMoveStart'')');
set(hMoveBut,'UserData',[],'ForegroundColor','k')
set(gcf,'WindowButtonMotionFcn',' ')
set(gcf,'WindowButtonUpFcn',' ')
set(gcf,'Pointer','Arrow')

return



elseif strcmp(lower(Action),lower('GraphicsSize'))
%=======================================================================
% Change size of next object clicked, provided it's editable and in figure
% "normal" mouse button decreases size
% "extend" mouse button increases size
% "alt" mouse button cancels operation
% spm_figure('GraphicsSize')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Resize: LeftMouse=shrink, MiddleMouse=grow, RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
set(F,'Pointer','Arrow')
set(F,'Name',tmp)

if strcmp(get(gcf,'SelectionType'),'alt')
	%-Quit option
	set(hBut,'ForegroundColor','k')
	return
end

h = gco(F);
NoDel=[F; findobj(F,'Tag','NoDelete')];
if ( ~any(h==NoDel) & (gcf==F) )
	u = strcmp(get(F,'SelectionType'),'extend');
	c  = get(h,'Type');
	if ( strcmp(c,'patch') | strcmp(c,'line') | strcmp(c,'image') )
		h = get(h,'Parent'); c = get(h,'Type'); end
	if strcmp(c,'text')
		set(h,'Fontsize',(get(h,'FontSize')+2*(2*u-1))); end
	if strcmp(c,'axes') | strcmp(c,'uicontrol')
		if u; u = 1.24; else u = 1/1.24; end
		P = get(h,'Position');
		set(h,'Position',[P(1:2)-P(3:4)*(u-1)/2, P(3:4)*u])
	end
end

set(hBut,'ForegroundColor','k')

return



elseif strcmp(lower(Action),lower('GraphicsText'))
%=======================================================================
% Add text annotation to a figure
% spm_figure('GraphicsText')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Select starting position, edit text widget. RightMouse=cancel..');
set(F,'Pointer','BotL')
waitforbuttonpress;
set(F,'Pointer','Arrow')
set(F,'Name',tmp)

if strcmp(get(gcf,'SelectionType'),'alt')
	%-Quit option
	set(hBut,'ForegroundColor','k')
	return
end

set(F,'Units','Normalized')
CPt = get(F,'CurrentPoint');

%-Set up axes for text widget
axes('Position',[CPt, 0.1, 0.1],'Visible','off');

%-Set uicontrol object and Callback
set(F,'Units','Pixels')
P = get(F,'Position'); P  = P(3:4);		% Figure dimensions {pixels}
uicontrol(F,'Style','Edit',...
	'Position',[CPt(1:2).*P, (1-CPt(1))*P(1), 20],...
	'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Left',...
	'Callback','text(0,0,get(gco,''String'')); delete(gco)');

set(hBut,'ForegroundColor','k')

return


elseif strcmp(lower(Action),lower('GraphicsTextEdit'))
%=======================================================================
% Add text annotation to a figure
% spm_figure('GraphicsText')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Select text to edit. RightMouse=cancel..');
set(F,'Pointer','Circle')
waitforbuttonpress;
set(F,'Pointer','Arrow')
set(F,'Name',tmp)

if strcmp(get(gcf,'SelectionType'),'alt')
	%-Quit option
	set(hBut,'ForegroundColor','k')
	return
end

NoDel=[F; findobj(F,'Tag','NoDelete')];
if any(gco==NoDel) | gcf~=F | ~strcmp(get(gco,'Type'),'text')
	set(hBut,'ForegroundColor','k')
	return
end

%-Save units of various objects
Funits = get(F,'Units');
Aunits = get(gca,'Units');
Tunits = get(gca,'Units');

%-Get locations
set(F,'Units','Pixels')
set(gca,'Units','Pixels')
set(gco,'Units','Pixels')
tmp = [1,1,0,0].*get(gca,'Position') ...
	+ [1,1,0,0].*get(gco,'Extent')...
	+ [0,0,max([0,0,1,1].*get(gco,'Extent')),20];

%-Create editable text widget to adjust text string
h = uicontrol(F,'Style','Edit',...
	'String',get(gco,'String'),...
	'Position',tmp,...
	'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Left',...
	'UserData',gco,...
	'Callback',['set(get(gco,''UserData''),'...
		'''String'',get(gco,''String'')), delete(gco)']);

%-Reset stuff
set(hBut,'ForegroundColor','k')
set(F,'Units',Funits)
set(gca,'Units',Aunits)
set(gco,'Units',Tunits)

return



else
%=======================================================================
error('Illegal Action string')


%=======================================================================
end
