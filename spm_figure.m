function varargout=spm_figure(varargin)
% Setup and callback functions for Graphics window
% FORMAT varargout=spm_figure(varargin)
%       - An embedded callback, multi-function function
%       - For detailed programmers comments, see format specifications
%         in main body of code
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
% objects or adding text). (This menu is also provided as a figure
% background "ContextMenu" - right-clicking on the figure background
% should bring up the menu.)
%
% Print: Creates a footnote (detailing the SPM version, user & date)
% and evaluates PRINTSTR (see spm_defaults.m). Graphics windows with
% multi-page axes are printed page by page.
%
% Clear: Clears the Graphics window. If in SPM usage (figure 'Tag'ed as
% 'Graphics') then all SPM windows are cleared and reset.
%
% Colormap options:
% * gray, hot, pink: Sets the colormap to its default values and loads
%                 either a grayscale, 'hot metal' or color map.
% * gray-hot, etc: Creates a 'split' colormap {128 x 3 matrix}.
%          The lower half is a gray scale and the upper half
%          is 'hot metal' or 'pink'.  This color map is used for
% 	   viewing 'rendered' SPMs on a PET, MRI or other background images
%
% Colormap effects:
% * Invert: Inverts (flips) the current color map.
% * Brighten and Darken: Brighten and Darken the current colourmap
% 	   using the MatLab BRIGHTEN command, with  beta's of +0.2 and -0.2
% 	   respectively.
%
% Editing: Right button ('alt' button) cancels operations
% * Cut  : Deletes the graphics object next selected (if deletable)
%          Select with middle mouse button to delete blocks of text,
%          or to delete individual elements from a plot.
% * Move : To re-position a text, uicontrol or axis object using a
%          'drag and drop' implementation (i.e. depress - move - release)
%          Using the middle 'extend' mouse button on a text object moves
%          the axes containing the text - i.e. blocks of text.
% * Size : Re-sizes the text, uicontrol or axis object next selected
%          {left button - decrease, middle button  - increase} by a factor
%          of 1.24 (or increase/decrease FontSize by 2 dpi)
% * Text : Creates an editable text widget that produces a text object as
%          its CallBack.
%          The text object is provided with a ContextMenu, obtained by
%          right-clicking ('alt') on the text, allowing text attributes
%          to be changed. Alternatively, the edit facilities on the window
%          menu bar or ContextMenu can be used.
% * Edit : To edit text, select a text object with the circle cursor,
%          and edit the text in the editable text widget that appears.
%          A middle 'extend' mouse click places a context menu on the text
%          object, facilitating easy modification of text atributes.
%
% For SPM usage, the figure should be 'Tag'ed as 'Graphics'.
%
% For SPM power users, and programmers, spm_figure provides utility
% routines for using the SPM graphics interface. Of particular use are
% the GetWin, FindWin and Clear functions See the embedded callback
% reference in the main body of spm_figure, below the help text.
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
% with a ToolBar and background context menu.
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
% FORMAT F = spm_figure('GetWin',Tag)
% Like spm_figure('FindWin',Tag), except that if 'Tag' is 'Graphics' or
% 'Interactive' and no such 'Tag'ged figure is found, one is created. Further,
% the "got" window is made current.
% Tag	- Figure 'Tag' to get, defaults to 'Graphics'
% F	- Figure number (if found/created) or empty (if not).
%
% FORMAT F = spm_figure('ParentFig',h)
% Finds window containing the object whose handle is specified
% h	- Handle of object whose parent figure is required
%	- If a vector, then first object handle is used
% F	- Number or parent figure
%
% FORMAT spm_figure('Clear',F,Tags)
% Clears figure, leaving ToolBar (& other objects with invisible handles)
% Optional third argument specifies 'Tag's of objects to delete.
% If figure F is 'Tag'ged 'Interactive' (SPM usage), then the window
% name and pointer are reset.
% F	- 'Tag' string or figure number of figure to clear, defaults to gcf
% Tags  - 'Tag's (string matrix or cell array of strings) of objects to delete
%         *regardless* of 'HandleVisibility'. Only these objects are deleted.
%         '!all' denotes all objects
%
% FORMAT str = spm_figure('DefPrintCmd')
% Returns default print command for SPM, as a string
%
% FORMAT spm_figure('Print',F)
% SPM print function: Appends footnote & executes PRINTSTR
% F	- [Optional] Figure to print. ('Tag' or figure number)
%	  Defaults to figure 'Tag'ed as 'Graphics'.
%	  If none found, uses CurrentFigure if avaliable.
% PRINTSTR - global variable holding print command to be evaluated
%	  Defaults to 'print -dps2 fig.ps'
% If objects 'Tag'ed 'NextPage' and 'PrevPage' are found, then the
% pages are shown and printed in order. In breif, pages are held as
% seperate axes, with ony one 'Visible' at any one time. The handles of
% the "page" axes are stored in the 'UserData' of the 'NextPage'
% object, while the 'PrevPage' object holds the current page number.
% See spm_help('!Disp') for details on setting up paging axes.
%
% FORMAT [hNextPage, hPrevPage, hPageNo] = spm_figure('NewPage',hPage)
% SPM pagination function: Makes objects with handles hPage paginated
% Creates pagination buttons if necessary.
% hPage	                        - Handles of objects to stick to this page
% hNextPage, hPrevPage, hPageNo - Handles of pagination controls
%
% FORMAT spm_figure('TurnPage',move,F)
% SPM pagination function: Turn to specified page
%
% FORMAT spm_figure('DeletePageControls',F)
% SPM pagination function: Deletes page controls
% F	- [Optional] Figure in which to attempt to turn the page
%         Defaults to 'Graphics' 'Tag'ged window
%
% FORMAT n = spm_figure('#page')
% Returns the current page number.
%
% FORMAT spm_figure('WaterMark',F,str,Tag,Angle,Perm)
% Adds watermark to figure windows.
% F	- Figure for watermark. Defaults to gcf
% str   - Watermark string. Defaults (missing or empty) to SPM
% Tag   - Tag for watermark axes. Defaults to ''
% Angle - Angle for watermark. Defaults to -45
% Perm  - If specified, then watermark is permanent (HandleVisibility 'off')
%
% FORMAT F = spm_figure('CreateWin',Tag,Name,Visible)
% Creates a full length WhiteBg figure 'Tag'ged Tag (if specified).
% F	  - Figure created
% Tag	  - Tag for window
% Name    - Name for window
% Visible - 'on' or 'off'
%
% FORMAT WS = spm_figure('GetWinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. WS=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
% (Function duplicated in spm.m, repeated to reduce inter-dependencies.)
%
% FORMAT FS = spm_figure('FontSizes',FS)
% Returns fontsizes FS scaled for the current display.
% FS     - (vector of) Font sizes to scale
%          [default [08,09,11,13,14,6:36]]
%
% FORMAT spm_figure('CreateBar',F)
% Creates toolbar in figure F (defaults to gcf). F can be a 'Tag'
% If the figure is 'Tag'ed as 'Graphics' (SPM usage), then the Print button
% callback is set to attempt to clear an 'Interactive' figure too.
%
% FORMAT h = spm_figure('FigContextMenu',F)
% Creates figure UIcontextMenu, with functionality of Figure ToolBar.
% F - handle of Figure for UIcontextMenu [default gcf]
% h - handle of UIContextMenu
%
% FORMAT h = spm_figure('TxtContextMenu',t)
% Creates UIcontextMenu for text objects
% t - handle of text object to ContextMenu, or figure to work in [default gcf]
% h - handle of UIcontextMenu
%
% FORMAT spm_figure('ColorMap')
% Callback for "ColorMap" buttons
%
% FORMAT h = spm_figure('GraphicsHandle',F)
% GUI choose object for handle identification. LeftMouse 'normal' returns
% handle, MiddleMouse 'extend' returns parents handle, RightMouse 'alt' cancels.
% F - figure to do a GUI "handle ID" in [Default gcbf]
%
% FORMAT spm_figure('GraphicsCut',F)
% Graphics cut - Callback for "Cut" button
% F - figure to do a GUI "cut" in [Default gcbf]
%
% FORMAT spm_figure('GraphicsMove',F)
% Graphics move - Callback for "Move" button
% F - figure to do a GUI "move" in [Default gcbf]
%
% FORMAT spm_figure('GraphicsMoveMotion')
% Callback for move function
%
% FORMAT spm_figure('GraphicsMoveEnd')
% Callback for move function
%
% FORMAT spm_figure('GraphicsReSize',F)
% Graphics resize - Callback for "Size" button
% F - figure to do a GUI "size" in [Default gcbf]
%
% FORMAT spm_figure('GraphicsText',F)
% Graphcis create text - Callback for "Text" button
% F - figure to do a GUI "text" in [Default gcbf]
%
% FORMAT spm_figure('GraphicsTextEdit',F)
% Graphics text edit - Callback for "Edit" button
% F - figure to do a GUI "edit" in [Default gcbf]
%_______________________________________________________________________


%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'Create'; else, Action = varargin{1}; end

switch lower(Action), case 'create'
%=======================================================================
% F = spm_figure('Create',Tag,Name,Visible)
%-Condition arguments
if nargin<4, Visible='on'; else, Visible=varargin{4}; end
if nargin<3, Name=''; else, Name=varargin{3}; end
if nargin<2, Tag=''; else, Tag=varargin{2}; end

F = spm_figure('CreateWin',Tag,Name,Visible);
spm_figure('CreateBar',F)
varargout = {F};


case 'findwin'
%=======================================================================
% F=spm_figure('FindWin',F)
% F=spm_figure('FindWin',Tag)
%-Find window: Find window with FigureNumber# / 'Tag' attribute
%-Returns empty if window cannot be found - deletes multiple tagged figs.

if nargin<2, F='Graphics'; else, F=varargin{2}; end

if isempty(F)
	% Leave F empty
elseif ischar(F)
	% Finds Graphics window with 'Tag' string - delete multiples
	Tag=F;
	F = findobj(get(0,'Children'),'Flat','Tag',Tag);
	if length(F) > 1
		% Multiple Graphics windows - close all but most recent
		close(F(2:end))
		F = F(1);
	end
else
	% F is supposed to be a figure number - check it
	if ~any(F==get(0,'Children')), F=[]; end
end
varargout = {F};

case 'getwin'
%=======================================================================
% F=spm_figure('GetWin',Tag)

if nargin<2, Tag='Graphics'; else, Tag=varargin{2}; end
F = spm_figure('FindWin',Tag);

if isempty(F)
	if ischar(Tag)
		switch Tag, case 'Graphics'
			F = spm_figure('Create','Graphics','Graphics');
		case 'Interactive'
			F = spm('CreateIntWin');
		end
	end
else
	set(0,'CurrentFigure',F);
end
varargout = {F};

case 'parentfig'
%=======================================================================
% F=spm_figure('ParentFig',h)
if nargin<2, error('No object specified'), else, h=varargin{2}; end
F = get(h(1),'Parent');
while ~strcmp(get(F,'Type'),'figure'), F=get(F,'Parent'); end
varargout = {F};


case 'clear'
%=======================================================================
% spm_figure('Clear',F,Tags)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<3, Tags=[]; else, Tags=varargin{3}; end
if nargin<2, F=get(0,'CurrentFigure'); else, F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Clear figure
%-----------------------------------------------------------------------
if isempty(Tags)
	%-Clear figure of objects with 'HandleVisibility' 'on'
	delete(findobj(get(F,'Children'),'flat','HandleVisibility','on'));
	%-Reset figures callback functions
	set(F,'KeyPressFcn','',...
		'WindowButtonDownFcn','',...
		'WindowButtonMotionFcn','',...
		'WindowButtonUpFcn','')
	%-If this is the 'Interactive' window, reset name & UserData
	if strcmp(get(F,'Tag'),'Interactive')
		set(F,'Name','','UserData',[]), end
else
	%-Clear specified objects from figure
	cSHH = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	if ischar(Tags); Tags=cellstr(Tags); end
	if any(strcmp(Tags(:),'!all'))
		delete(get(F,'Children'))
	else
	    for tag = Tags(:)'
		delete(findobj(get(F,'Children'),'flat','Tag',tag{:}));
	    end
	end	
	set(0,'ShowHiddenHandles',cSHH)
end
set(F,'Pointer','Arrow')


case 'defprintcmd'
%=======================================================================
% spm_figure('DefPrintCmd')
varargout = {'print -dpsc2 -painters -append -noui '};


case 'print'
%=======================================================================
% spm_figure('Print',F)

%-Arguments & defaults
if nargin<2, F='Graphics'; else, F=varargin{2}; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
F=spm_figure('FindWin',F);
if isempty(F), F = get(0,'CurrentFigure'); end
if isempty(F), return, end

%-Note current figure, & switch to figure to print
cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',F)

%-See if window has paging controls
hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
iPaged    = ~isempty(hNextPage);

%-Construct print command
%-----------------------------------------------------------------------
if any(strcmp(who('global'),'PRINTSTR'))
	global PRINTSTR
else
	PRINTSTR = [spm_figure('DefPrintCmd'),'spm.ps'];
end

%-Create footnote with SPM version, username, date and time.
%-----------------------------------------------------------------------
FNote = sprintf('%s%s: %s',spm('ver'),spm('GetUser',' (%s)'),spm('time'));
%-Delete old tag lines, and print new one
delete(findobj(F,'Tag','SPMprintFootnote'));
axes('Position',[0.005,0.005,0.1,0.1],...
	'Visible','off',...
	'Tag','SPMprintFootnote')
text(0,0,FNote,'FontSize',6);


%-Temporarily change all units to normalized prior to printing
% (Fixes bizzarre problem with stuff jumping around!)
%-----------------------------------------------------------------------
H  = findobj(get(F,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
set(H,'Units','normalized')


%-Print
%-----------------------------------------------------------------------
err = 0;
if ~iPaged
	try, eval(PRINTSTR), catch, err=1; end
else
	hPg       = get(hNextPage,'UserData');
	Cpage     = get(hPageNo,  'UserData');
	nPages    = size(hPg,1);

	set([hNextPage,hPrevPage,hPageNo],'Visible','off')
	if Cpage~=1
		set(hPg{Cpage,1},'Visible','off'), end
	for p = 1:nPages
		set(hPg{p,1},'Visible','on')
		try, eval(PRINTSTR), catch, err=1; end
		set(hPg{p,1},'Visible','off')
	end
	set(hPg{Cpage,1},'Visible','on')
	set([hNextPage,hPrevPage,hPageNo],'Visible','on')
end

if err
	errstr = lasterr;
	tmp = [find(abs(errstr)==10),length(errstr)+1];
	str = {errstr(1:tmp(1)-1)};
	for i = 1:length(tmp)-1
		if tmp(i)+1 < tmp(i+1) 
			str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
		end
	end
	str = {str{:},	'','- print command is:',['    ',PRINTSTR],...
			'','- current directory is:',['    ',pwd],...
			'','            * nothing has been printed *'};
	spm('alert!',str,'printing problem...',sqrt(-1));
end

set(H,{'Units'},un)
set(0,'CurrentFigure',cF)


case 'newpage'
%=======================================================================
% [hNextPage, hPrevPage, hPageNo] = spm_figure('NewPage',h)
if nargin<2 | isempty(varargin{2}), error('No handles to paginate')
	else, h=varargin{2}(:)'; end

%-Work out which figure we're in
F = spm_figure('ParentFig',h(1));

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

%-Create pagination widgets if required
%-----------------------------------------------------------------------
if isempty(hNextPage)
	WS = spm('WinScale');
	FS = spm('FontSizes');
	hNextPage = uicontrol(F,'Style','Pushbutton',...
		'HandleVisibility','on',...
		'String','>','FontSize',FS(10),...
		'ToolTipString','next page',...
		'Callback','spm_figure(''TurnPage'',''+1'',gcbf)',...
		'Position',[580 020 015 015].*WS,...
		'ForegroundColor',[0 0 0],...
		'Tag','NextPage','UserData',[]);
	hPrevPage = uicontrol(F,'Style','Pushbutton',...
		'HandleVisibility','on',...
		'String','<','FontSize',FS(10),...
		'ToolTipString','previous page',...
		'Callback','spm_figure(''TurnPage'',''-1'',gcbf)',...
		'Position',[565 020 015 015].*WS,...
		'Visible','on',...
		'Enable','off',...
		'Tag','PrevPage');
	hPageNo = uicontrol(F,'Style','Text',...
		'HandleVisibility','on',...
		'String','1',...
		'FontSize',FS(6),...
		'HorizontalAlignment','center',...
		'BackgroundColor','w',...
		'Position',[550 005 060 015].*WS,...
		'Visible','on',...
		'UserData',1,...
		'Tag','PageNo','UserData',1);
end

%-Add handles for this page to UserData of hNextPage
%-Make handles for this page invisible if PageNo>1
%-----------------------------------------------------------------------
mVis    = strcmp('on',get(h,'Visible'));
hPg     = get(hNextPage,'UserData');
if isempty(hPg)
	hPg = {h(mVis), h(~mVis)};
else
	hPg = [hPg; {h(mVis), h(~mVis)}];
	set(h(mVis),'Visible','off')
end
set(hNextPage,'UserData',hPg)

%-Return handles to pagination controls if requested
if nargout>0, varargout = {[hNextPage, hPrevPage, hPageNo]}; end


case 'turnpage'
%=======================================================================
% spm_figure('TurnPage',move,F)
if nargin<3, F='Graphics'; else, F=varargin{3}; end
if nargin<2, move=1; else, move=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
if isempty(hNextPage), return, end
hPg       = get(hNextPage,'UserData');
Cpage     = get(hPageNo,  'UserData');
nPages    = size(hPg,1);

%-Sort out new page number
if ischar(move), Npage = Cpage+eval(move); else, Npage = move; end
Npage = max(min(Npage,nPages),1);

%-Make current page invisible, new page visible, set page number string
set(hPg{Cpage,1},'Visible','off')
set(hPg{Npage,1},'Visible','on')
set(hPageNo,'UserData',Npage,'String',sprintf('%d / %d',Npage,nPages))

%-Disable appropriate page turning control if on first/last page (for neatness)
if Npage==1, set(hPrevPage,'Enable','off')
else, set(hPrevPage,'Enable','on'), end
if Npage==nPages, set(hNextPage,'Enable','off')
else, set(hNextPage,'Enable','on'), end



case 'deletepagecontrols'
%=======================================================================
% spm_figure('DeletePageControls',F)
if nargin<2, F='Graphics'; else, F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

delete([hNextPage hPrevPage hPageNo])


case '#page'
%=======================================================================
% n = spm_figure('#Page',F)
if nargin<2, F='Graphics'; else, F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
if isempty(hNextPage)
	n = 1;
else
	n = size(get(hNextPage,'UserData'),1)+1;
end
varargout = {n};


case 'watermark'
%=======================================================================
% spm_figure('WaterMark',F,str,Tag,Angle,Perm)
if nargin<6, HVis='on'; else, HVis='off'; end
if nargin<5, Angle=-45; else, Angle=varargin{5}; end
if nargin<4 | isempty(varargin{4}), Tag = 'WaterMark'; else, Tag=varargin{4}; end
if nargin<3 | isempty(varargin{3}), str = 'SPM';       else, str=varargin{3}; end
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Specify watermark color from background colour
%-----------------------------------------------------------------------
Colour = get(F,'Color');
%-Only mess with grayscale backgrounds
if ~all(Colour==Colour(1)), return, end
%-Work out colour - lighter unless grey value > 0.9
Colour = Colour+(2*(Colour(1)<0.9)-1)*0.02;

cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',F)
Units=get(F,'Units');
set(F,'Units','normalized');
h = axes('Position',[0.45,0.5,0.1,0.1],...
	'Units','normalized',...
	'Visible','off',...
	'Tag',Tag);
set(F,'Units',Units)
text(0.5,0.5,str,...
	'FontSize',spm('FontSize',80),...
	'FontWeight','Bold',...
	'FontName',spm_platform('Font','times'),...
	'Rotation',Angle,...
	'HorizontalAlignment','Center',...
	'VerticalAlignment','middle',...
	'EraseMode','normal',...
	'Color',Colour,...
	'ButtonDownFcn',[...
		'if strcmp(get(gcbf,''SelectionType''),''open''),',...
			'delete(get(gcbo,''Parent'')),',...
		'end'])
set(h,'HandleVisibility',HVis)
set(0,'CurrentFigure',cF)


case 'createwin'
%=======================================================================
% F=spm_figure('CreateWin',Tag,Name,Visible)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<4 | isempty(varargin{4}), Visible='on'; else, Visible=varargin{4}; end
if nargin<3, Name=''; else, Name = varargin{3}; end
if nargin<2, Tag='';  else, Tag  = varargin{2}; end

WS   = spm('WinScale');				%-Window scaling factors
FS   = spm('FontSizes');			%-Scaled font sizes
PF   = spm_platform('fonts');			%-Font names (for this platform)
Rect = spm('WinSize','Graphics','raw').*WS;	%-Graphics window rectangle

F      = figure(...
	'Tag',Tag,...
	'Position',Rect,...
	'Resize','off',...
	'MenuBar','none',...
	'Color','w',...
	'ColorMap',gray(64),...
	'DefaultTextColor','k',...
	'DefaultTextInterpreter','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
	'DefaultAxesColor','w',...
	'DefaultAxesXColor','k',...
	'DefaultAxesYColor','k',...
	'DefaultAxesZColor','k',...
	'DefaultAxesFontName',PF.helvetica,...
	'DefaultPatchFaceColor','k',...
	'DefaultPatchEdgeColor','k',...
	'DefaultSurfaceEdgeColor','k',...
	'DefaultLineColor','k',...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(10),...
	'DefaultUicontrolInterruptible','on',...
	'PaperType','A4',...
	'PaperUnits','normalized',...
	'PaperPosition',[.0726 .0644 .854 .870],...
	'InvertHardcopy','off',...
	'Renderer','zbuffer',...
	'Visible','off');
if ~isempty(Name)
	set(F,'Name',sprintf('%s%s: %s',spm('ver'),...
		spm('GetUser',' (%s)'),Name),'NumberTitle','off')
end
spm_figure('FigContextMenu',F);
set(F,'Visible',Visible)
varargout = {F};


case 'getwinscale'
%=======================================================================
% WS = spm_figure('GetWinScale')
warning('spm_figure(''GetWinScale''... is Grandfathered: use spm(''WinScale''')
varargout = {spm('WinScale')};


case 'fontsizes'
%=======================================================================
% FS = spm_figure('FontSizes',FS)
warning('spm_figure(''FontSizes''... is Grandfathered: use spm(''FontSizes''')
if nargin<2, FS=[08,09,11,13,14,6:36]; else, FS=varargin{2}; end
varargout = {round(FS*min(spm('WinScale')))};


case 'createbar'
%=======================================================================
% spm_figure('CreateBar',F)
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Get position and size parameters
%-----------------------------------------------------------------------
cUnits = get(F,'Units');
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
FS    = round(10*min(S_Gra));			% uicontrol font size

%-Delete any existing 'ToolBar' 'Tag'ged objects
%-----------------------------------------------------------------------
cSHH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
delete(findobj(F,'Tag','ToolBar'));
set(0,'ShowHiddenHandles',cSHH)

%-Create Frame for controls
%-----------------------------------------------------------------------
uicontrol(F,'Style', 'Frame',...
	'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
	'Tag','ToolBar','HandleVisibility','callback');

%-Create uicontrol objects
%-----------------------------------------------------------------------
uicontrol(F,'String','Print','ToolTipString','print figure',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''Print'',gcf)',...
	'FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','b'); x = x+sx;

h = uicontrol(F,'String','Clear','ToolTipString','clear figure',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''Clear'',gcbf)',...
	'FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
        'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','b'); x = x+sx+dx;
if strcmp(get(F,'Tag'),'Graphics')	%-Do a full SPM interface clear
	set(h,'CallBack','spm(''Clear'',''Interactive'',gcbf)',...
		'ToolTipString','clear figure & reset SPM GUI')
end

uicontrol(F,'Style','PopUp','String',...
	'ColorMap|gray|hot|pink|gray-hot|gray-pink',...
	'ToolTipString','change colormap',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'spm_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'FontSize',FS,...
	'Tag','ToolBar','HandleVisibility','callback',...
	'UserData','Pop'); x = x + 2*sx;

uicontrol(F,'Style','PopUp','String','Effects|invert|brighten|darken',...
	'ToolTipString','colormap effects',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'spm_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'FontSize',FS,...
	'Tag','ToolBar','HandleVisibility','callback',...
	'UserData','Pop'); x = x + 2*sx + dx;

uicontrol(F,'String','cut','ToolTipString','delete a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsCut'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','move','ToolTipString','move a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsMove'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','resize','ToolTipString','resize a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsReSize'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','text','ToolTipString','create text annotation',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsText'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','edit','ToolTipString','edit a text string',...
	'Position',[x y sx sy],...
	'CallBack','spm_figure(''GraphicsTextEdit'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx+dx;

uicontrol(F,'String','?','ToolTipString','spm_figure help',...
	'Position',[x y dx sy],...
	'CallBack','spm_help(''spm_figure.m'')','FontSize',FS,...
	'Interruptible','off','BusyAction','queue',...
	'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','g'); x = x+2*dx;

set(F,'Units',cUnits)


case 'figcontextmenu'
%=======================================================================
% h = spm_figure('FigContextMenu',F)
if nargin<2
	F = get(0,'CurrentFigure');
	if isempty(F), error('no figure'), end
else
	F = spm_figure('FindWin',varargin{2});
	if isempty(F), error('no such figure'), end
end

h = uicontextmenu('Parent',F,'HandleVisibility','CallBack');
uimenu(h,'Label','Print','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''Print'',gcf)')
uimenu(h,'Label','Clear','HandleVisibility','CallBack',...
	    'CallBack','spm_figure(''Clear'',gcbf)');
%uimenu(h,'Label','Clear','HandleVisibility','CallBack',...
%	'CallBack','spm(''Clear'',''Interactive'',gcbf)');

hC = uimenu(h,'Label','Colormap','HandleVisibility','CallBack','Separator','on');
uimenu(hC,'Label','gray','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''gray'')')
uimenu(hC,'Label','hot','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''hot'')')
uimenu(hC,'Label','pink','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''pink'')')
uimenu(hC,'Label','gray-hot','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''gray-hot'')')
uimenu(hC,'Label','gray-pink','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''gray-pink'')')

hE = uimenu(h,'Label','Effects','HandleVisibility','CallBack');
uimenu(hE,'Label','invert','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''invert'')')
uimenu(hE,'Label','brighten','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''brighten'')')
uimenu(hE,'Label','darken','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''ColorMap'',''darken'')')

uimenu(h,'Label','cut','HandleVisibility','CallBack','Separator','on',...
	'CallBack','spm_figure(''GraphicsCut'')')
uimenu(h,'Label','move','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''GraphicsMove'')')
uimenu(h,'Label','resize','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''GraphicsReSize'')')
uimenu(h,'Label','text','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''GraphicsText'')')
uimenu(h,'Label','edit','HandleVisibility','CallBack',...
	'CallBack','spm_figure(''GraphicsTextEdit'')')

uimenu(h,'Label','help','HandleVisibility','CallBack','Separator','on',...
	'CallBack','spm_help(''spm_figure.m'')')

%uimenu(h,'Label','Identify handle','HandleVisibility','CallBack','Separator','on',...
%	'CallBack','spm_figure(''GraphicsHandle'')')

set(F,'UIContextMenu',h)
varargout = {h};


case 'txtcontextmenu'
%=======================================================================
% h = spm_figure('TxtContextMenu',t)
if nargin<2
	F = get(0,'CurrentFigure'); t=[];
	if isempty(F), error('no figure'), end
elseif ischar(varargin{2})		%-Figure 'Tag'
	F = spm_figure('FindWin',varargin{2}); t=[];
	if isempty(F), error(sprintf('no figure ''Tag''ged ''%s''',...
		varargin{2})), end
elseif ishandle(varargin{2})	%-handle: figure or text object
	switch get(varargin{2},'Type')
	case 'figure'
		F = varargin{2};
		t = [];
	case 'text'
		t = varargin{2};
		F = spm_figure('ParentFig',t);
	otherwise
		error('handle not a figure or text object')
	end
else
	error('invalid handle')
end



h = uicontextmenu('Parent',F);
uimenu(h,'Label','Delete','CallBack','delete(gco)')
uimenu(h,'Label','Move','CallBack','spm_figure(''GraphicsMove'')')

tmp = uimenu(h,'Label','Font','Separator','on');
uimenu(tmp,	'Label','normal','CallBack','set(gco,''FontAngle'',''normal'')')
uimenu(tmp,	'Label','italic','CallBack','set(gco,''FontAngle'',''italic'')')
uimenu(tmp,	'Label','oblique',...
		'CallBack','set(gco,''FontAngle'',''oblique'')')

uimenu(tmp,'Separator','on',...
		'Label','Helvetica','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''helvetica''))')
uimenu(tmp,	'Label','Times','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''times''))')
uimenu(tmp,	'Label','Courier','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''courier''))')
uimenu(tmp,	'Label','Symbol','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''symbol''))')

uimenu(tmp,'Separator','on',...
		'Label','light','CallBack','set(gco,''FontWeight'',''light'')')
uimenu(tmp,	'Label','normal','CallBack','set(gco,''FontWeight'',''normal'')')
uimenu(tmp,	'Label','demi','CallBack','set(gco,''FontWeight'',''demi'')')
uimenu(tmp,	'Label','bold','CallBack','set(gco,''FontWeight'',''bold'')')

tmp = uimenu(h,'Label','FontSize');
uimenu(tmp,	'Label', '6','CallBack','set(gco,''FontSize'', 6)')
uimenu(tmp,	'Label', '8','CallBack','set(gco,''FontSize'', 8)')
uimenu(tmp,	'Label','10','CallBack','set(gco,''FontSize'',10)')
uimenu(tmp,	'Label','12','CallBack','set(gco,''FontSize'',12)')
uimenu(tmp,	'Label','14','CallBack','set(gco,''FontSize'',14)')
uimenu(tmp,	'Label','16','CallBack','set(gco,''FontSize'',16)')
uimenu(tmp,	'Label','18','CallBack','set(gco,''FontSize'',18)')
uimenu(tmp,	'Label','20','CallBack','set(gco,''FontSize'',20)')
uimenu(tmp,	'Label','24','CallBack','set(gco,''FontSize'',24)')
uimenu(tmp,	'Label','28','CallBack','set(gco,''FontSize'',28)')
uimenu(tmp,	'Label','36','CallBack','set(gco,''FontSize'',36)')
uimenu(tmp,	'Label','48','CallBack','set(gco,''FontSize'',48)')
uimenu(tmp,	'Label','64','CallBack','set(gco,''FontSize'',64)')
uimenu(tmp,	'Label','96','CallBack','set(gco,''FontSize'',96)')

tmp = uimenu(h,'Label','Rotation');
uimenu(tmp,	'Label',  '0','CallBack','set(gco,''Rotation'',0)')
uimenu(tmp,	'Label', '90','CallBack','set(gco,''Rotation'',90)')
uimenu(tmp,	'Label','180','CallBack','set(gco,''Rotation'',180)')
uimenu(tmp,	'Label','-90','CallBack','set(gco,''Rotation'',270)')

uimenu(h,	'Label','Edit','CallBack','set(gco,''Editing'',''on'')')

tmp = uimenu(h,	'Label','Interpreter');
uimenu(tmp,	'Label','none','CallBack','set(gco,''Interpreter'',''none'')')
uimenu(tmp,	'Label','TeX','CallBack','set(gco,''Interpreter'',''tex'')')

uimenu(h,'Separator','on','Label','Get handle','CallBack','gco')

if ~isempty(t), set(t,'UIContextMenu',h), end
varargout = {h};


case 'colormap'
%=======================================================================
% spm_figure('Colormap',ColAction,h)
if nargin<3, h=[]; else, h=varargin{3}; end
if nargin<2, ColAction='gray'; else, ColAction=varargin{2}; end
if ~ischar(ColAction)
	if ColAction==1, return, end
	Actions   = get(gcbo,'String');
	ColAction = deblank(Actions(ColAction,:));
end

switch lower(ColAction), case 'gray'
	colormap(gray(64))
case 'hot'
	colormap(hot(64))
case 'pink'
	colormap(pink(64))
case 'gray-hot'
	tmp = hot(64 + 16);  tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'gray-pink'
	tmp = pink(64 + 16); tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'invert'
	colormap(flipud(colormap))
case 'brighten'
	colormap(brighten(colormap, 0.2))
case 'darken'
	colormap(brighten(colormap, -0.2))
otherwise
	error('Illegal ColAction specification')
end

case 'graphicshandle'
%=======================================================================
% h = spm_figure('GraphicsHandle',F)
if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

tmp = get(F,'Name');
set(F,'Name',...
	'Handle: Select item to identify, MiddleMouse=parent, RightMouse=cancel...');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(SelnType,'alt') & ~isempty(h) & gcf==F
	str = sprintf('Selected (%s) object',get(h,'Type'));
	if strcmp(SelnType,'normal')
		str = sprintf('%s: handle',str);
	else
		h = get(h,'Parent');
		str = sprintf('%s: handle of parent (%s) object',str,get(h,'Type'));
	end
	if nargout==0
		assignin('base','ans',h)
		fprintf('\n%s: \n',str)
		ans = h
	else
		varargout={h};
	end
else
	varargout={[]};
end



case 'graphicscut'
%=======================================================================
% spm_figure('GraphicsCut',F)
% Delete next object clicked, provided it's deletable, in this gcbf & not UIcon
% "normal" mouse button deletes object (or parent if image, line, patch, surface)
% "extend" mouse button deletes parent (if text), or object (if image, line,
%          patch, surface)
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name','Cut: Select object to delete. RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if strcmp(get(h,'HandleVisibility'),'on') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hType,{'root','figure','uimenu','uicontrol'}) & ...
		gcf==F
	if strcmp(hType,'axes')
		delete(h)
	elseif strcmp(hType,'text')
		if strcmp(SelnType,'extend'), delete(get(h,'Parent'))
			else, delete(h), end
	elseif any(strcmp(hType,{'image','line','patch','surface'}))
		if strcmp(SelnType,'extend'), delete(h)
			else, delete(get(h,'Parent')), end
	end
end
set(hBut,'ForegroundColor','k')


case 'graphicsmove'
%=======================================================================
% spm_figure('GraphicsMove',F)
% Move the next object clicked, provided it's movable and in this figure
% "normal" mouse button moves object (or parent if image, line, patch, surface)
% "extend" mouse button moves parent (if text)
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hMoveBut   = gcbo;
set(hMoveBut,'ForegroundColor','r')
tmp        = get(F,'Name');
set(F,'Name','Move: MiddleMouse moves parent. RightMouse=cancel');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
hPress     = gco(F);
hPressType = get(hPress,'Type');
SelnType   = get(gcf,'SelectionType');
set(F,'Pointer','Fleur','Name',tmp)

if ~strcmp(get(hPress,'HandleVisibility'),'off') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hPressType,{'root','figure','uimenu','uicontrol'})&...
		gcf==F
	MS.cFUnits = get(F,'Units');
	set(F,'Units','Pixels')
	MS.OPt  = get(F,'CurrentPoint');

	if ( strcmp(SelnType,'extend') & strcmp(hPressType,'text') ) | ...
			any(strcmp(hPressType,{'image','line','patch','surface'}))
		hMove = get(hPress,'Parent');
	elseif any(strcmp(hPressType,{'axes','text'}))
		hMove = hPress;
	else
		set(hMoveBut,'ForegroundColor','k')
		set(F,'Pointer','Arrow')
		return
	end

	%-Store info in UserData of hPress
	MS.hMove    = hMove;
	MS.hMoveBut = hMoveBut;
	MS.chMUnits = get(hMove,'Units');
	set(hMove,'Units','Pixels');
	MS.OPos     = get(hMove,'Position');
	MS.UserData = get(F,'UserData');
	set(hPress,'UserData',MS)

	%-Set Motion & ButtonUp functions
	set(F,'WindowButtonMotionFcn','spm_figure(''GraphicsMoveMotion'')')
	set(F,'WindowButtonUpFcn','spm_figure(''GraphicsMoveEnd'')')
else
	set(hMoveBut,'ForegroundColor','k')
	set(F,'Pointer','Arrow')
end


case 'graphicsmovemotion'
%=======================================================================
% spm_figure('GraphicsMoveMotion')
MS = get(gco,'UserData');
set(MS.hMove,'Units','Pixels',...
	'Position',...
	MS.OPos + [get(gcf,'CurrentPoint')-MS.OPt(1:2),0*MS.OPos(3:end)])


case 'graphicsmoveend'
%=======================================================================
% spm_figure('GraphicsMoveEnd')
MS = get(gco,'UserData');
hType = get(gco,'Type');
if any(strcmp(hType,{'image','line','patch','surface'}))
	set(get(gco,'Parent'),'Units',MS.chMUnits);
	set(gco,'UserData',MS.UserData);
else
	set(gco,'Units',MS.chMUnits,...
		'UserData',MS.UserData);
end

set(gcf,'Units',MS.cFUnits,...
	'WindowButtonMotionFcn','',...
	'WindowButtonUpFcn','',...
	'Pointer','Arrow')
set(MS.hMoveBut,'ForegroundColor','k')


case 'graphicsresize'
%=======================================================================
% spm_figure('GraphicsReSize',F)
% Change size of next object clicked, provided it's editable and in figure
% "normal" mouse button decreases size
% "extend" mouse button increases size
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp   = get(F,'Name');
set(F,'Name',...
	'Resize: LeftMouse=shrink, MiddleMouse=grow, RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
u = 2*strcmp(SelnType,'extend')-1;
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(get(h,'HandleVisibility'),'off') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hType,{'root','figure','uimenu','uicontrol'}) & ...
		gcf==F
	if any(strcmp(hType,{'image','line','patch','surface'}))
		h = get(h,'Parent'); hType = get(h,'Type'); end
	if strcmp(hType,'text')
		set(h,'Fontsize',(get(h,'FontSize')+2*u))
	elseif strcmp(hType,'axes')
		if u==1; u = 1.24; else, u = 1/1.24; end
		P = get(h,'Position');
		set(h,'Position',[P(1:2)-P(3:4)*(u-1)/2, P(3:4)*u])
	end
end
set(hBut,'ForegroundColor','k')



case 'graphicstext'
%=======================================================================
% spm_figure('GraphicsText',F)
% Add text annotation to a figure

if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Select starting position, edit text widget. RightMouse=cancel..');
set(F,'Pointer','BotL')
waitforbuttonpress;
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(get(gcf,'SelectionType'),'alt') & gcf==F
	cUnits = get(F,'Units');
	set(F,'Units','Normalized')
	CPt = get(F,'CurrentPoint');
	
	%-Set up axes for text widget - draw empty text object
	axes('Position',[CPt, 0.1, 0.1],'Visible','off');
	h = text(0,0,'','UIContextMenu',spm_figure('TxtContextMenu'));
	
	%-Set uicontrol object and Callback
	set(F,'Units','Pixels')
	P = get(F,'Position');
	uicontrol(F,'Style','Edit',...
		'ToolTipString','enter text here',...
		'Position',[CPt(1:2).*P(3:4), (1-CPt(1))*P(3), 22],...
		'BackGroundColor',[0.8,0.8,1.0],...
		'HorizontalAlignment','Left',...
		'UserData',h,...
		'Callback',['set(get(gcbo,''UserData''),'...
			'''String'',get(gcbo,''String'')), delete(gcbo)']);
	set(F,'Units',cUnits)
end
set(hBut,'ForegroundColor','k')


case 'graphicstextedit'
%=======================================================================
% spm_figure('GraphicsTextEdit',F)
% Edit text annotation to a figure

if nargin<2, F=gcbf; else, F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut     = gcbo;
set(hBut,'ForegroundColor','r')
tmp      = get(F,'Name');
set(F,'Name','Select text to edit. RightMouse=cancel..');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if strcmp(get(h,'HandleVisibility'),'off') | ...
		strcmp(SelnType,'alt') | ...
		~strcmp(get(h,'Type'),'text') | ...
		gcf ~= F
	set(hBut,'ForegroundColor','k')
	return
end

if strcmp(SelnType,'extend')
	set(hBut,'ForegroundColor','k')
	set(h,'UIContextMenu',spm_figure('TxtCOntextMenu'))
	return
end

%-Save units of various objects
cFUnits     = get(F,'Units');
cAUnits     = get(gca,'Units');
chUnits     = get(h,'Units');
chFontUnits = get(h,'FontUnits');

%-Get locations
set(F,'Units','Pixels')
set(gca,'Units','Pixels')
set(h,'Units','Pixels')
set(h,'FontUnits','Points')
tExtent = get(h,'Extent');
tmp = [1,1,0,0].*get(gca,'Position') ...
	+ [1,1,0,0].*tExtent...
	+ [0,0,1.2*max([0,0,1,1].*tExtent),ceil(22*get(h,'FontSize')/12)];

%-Create editable text widget to adjust text string
uicontrol(F,'Style','Edit',...
	'ToolTipString','edit text & press return',...
	'String',get(h,'String'),...
	'FontAngle',get(h,'FontAngle'),...
	'FontName',get(h,'FontName'),...
	'FontSize',get(h,'FontSize'),...
	'Position',tmp,...
	'BackGroundColor',[0.8,0.8,1.0],...
	'HorizontalAlignment',get(h,'HorizontalAlignment'),...
	'UserData',h,...
	'Callback',['set(get(gcbo,''UserData''),'...
		'''String'',get(gcbo,''String'')), delete(gcbo)']);

%-Reset stuff
set(hBut,'ForegroundColor','k')
set(F,'Units',cFUnits)
set(gca,'Units',cAUnits)
set(h,'Units',chUnits)
set(h,'FontUnits',chFontUnits)

otherwise
%=======================================================================
warning(['Illegal Action string: ',Action])


%=======================================================================
end
