function R1=spm_figure(Action,P2,P3,P4)
% Setup and callback functions for Graphics window
% FORMAT R1=spm_figure(Action,P2,P3,P4)
%	- An embedded callback, multi-function function
%       - For detailed programmers comments, see the format specifications
%         below
%_______________________________________________________________________
%
% spm_figure creates and manages the 'Results' window, referred to as
% the Graphics window. This window and these facilities may be used
% independently of SPM, and any number of Graphics windows my be used
% within the same MatLab session. (Though only one SPM 'Graphics' 'Tag'ed
% window is permitted.
%
% The Graphics window is provided with a menu bar at the top that
% facilitates editing and printing of the current graphic display,
% faciliating interactive editing of graphic output prior to printing
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
% Changing The Colormap:
% * Gray and Hot: Sets the colormap to its default values and loads
%                 either a grayscale or 'hot metal' color map.
% * Split: Loads a 'split' color map from 'Split.mat' {128 x 3
%          matrix}.  The lower half is a gray scale and the upper half
%          is 'hot metal'.  This color map is used for viewing 'rendered' 
%          SPM{Z} on a PET, MRI or other background images
%
% Invert: Inverts (flips) the current color map.
%
% Editing:
% * Cut  : Deletes the graphics object next selected (if deletable)
% * Move : To re-position a text, uicontrol or axis object using a
%          'drag and drop' implementation (i.e. depress - move - release)
% * Size : Re-sizes the text, uicontrol or axis object next selected
%          {left button - decrease, right button  - increase} by a factor
%          of 1.24 (or increase/decrease FontSize by 2 dpi)
% * Text : Creates an editable text widget that produces a text object as
%          its CallBack.
%          This text object can then be manipulated using the edit facilities.
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
% below.
%
% See also: spm_print, spm_clf
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%
%
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
% FORMAT F = spm_figure('Create',Tag)
% Create a full length WhiteBg figure 'Tag'ed Tag (if specified),
% with a ToolBar.
% Equivalent to spm_figure('CreateWin','Tag') and spm_figure('CreateBar')
% Tag	- 'Tag' string for figure.
% F	- Figure used
%
% FORMAT F = spm_figure('FindWin',F)
% Finds window with 'Tag' or figure numnber F - returns empty F if not found
% F	- (Input)  Figure to use [Optional] - 'Tag' string or figure number.
%	- Defaults to 'Graphics'
% F	- (Output) Figure number (if found) or empty (if not).
%
% FORMAT spm_figure('Clear',F)
% F	- 'Tag' string or figure number of figure to clear, defaults to gcf
% Clears figure, leaving ToolBar (& other objects 'Tag'ed as 'NoDelete')
% intact. If figure F is 'Tag'ged 'Interactive' (SPM usage), then the window
% name and pointer are reset.
%
% FORMAT F = spm_figure('CreateWin',Tag,Visible)
% Creates a full length WhiteBg figure 'Tag'ged Tag (if specified).
% F	  - Figure created
% Tag	  - Tag for window
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
% See spm_help('Disp') for details on setting up paging axes.
%
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
%____________________________________________________________________________


%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'Create'; nargin=1; end
if ~isstr(Action), P2=Action; Action='CreateBar'; end

if strcmp(Action,'Create')
%=======================================================================
% F = spm_figure('Create',Tag)
%-Condition arguments
Tag = [];
if (nargin==2), if isstr(P2), Tag=P2; end, end

F = spm_figure('CreateWin',Tag);

spm_figure('CreateBar',F)
R1 = F;
return


elseif strcmp(Action,'CreateWin')
%=======================================================================
% F=spm_figure('CreateWin',Tag,Visible)
if (nargin<3), Visible='on'; else, Visible=P3; end
if (nargin<2), Tag=[]; else, Tag=P2; end
S0     = get(0,'ScreenSize');
Fsca   = [S0(3)/1152 S0(4)/900]; Fsca = [Fsca,Fsca];
S_Gra  = [515 008 600 865].*Fsca;
F      = figure('Name','Results',...
	'NumberTitle','off',...
	'Position',S_Gra,...
	'Resize','off',...
	'Visible',Visible,...
	'PaperPosition',[.75 1.5 7 9.5]);
if isstr(Tag), set(F,'Tag',Tag), end
whitebg(F,'w')
colormap gray
% set(F,'DefaultTextFontSize',2*round(12*min(Fsca)/2));
R1 = F;
return


elseif strcmp(Action,'FindWin')
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
	if ~any(F==findobj(get(0,'Children'))), F=[]; end
end

R1 = F;
return



elseif strcmp(Action,'Clear')
%=======================================================================
% spm_figure('Clear',F)
%-Clear window, leaving 'NoDelete' 'Tag'ed objects, reset pointer & name

%-Sort out arguments
%-----------------------------------------------------------------------
if (nargin<2)
	if any(get(0,'Children')), F = gcf; else, F=[], end
else
	F = P2;
end
F=spm_figure('FindWin',F);
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


elseif strcmp(Action,'Print')
%=======================================================================
% spm_figure('Print',F,PFile)

%-Arguments & defaults
if nargin<3, PFile='fig.ps'; else, PFile=P3; end
if nargin<2, F='Graphics'; else, F=P2; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
F=spm_figure('FindWin',F);
if isempty(F), if any(get(0,'Children')), F = gcf; end, end
if isempty(F), return, end

%-See if window has paging controls
hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
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
axes('Position',[0.01,0.01,1,1],...
	'Visible','off',...
	'Tag','SPMprintFootnote')
text(0,0,FNote,'FontSize',10);

%-Print
%-----------------------------------------------------------------------
if ~iPaged
	eval(PrintCmd)
else
	hAxes     = get(hNextPage,'UserData');
	Cpage     = get(hPrevPage,'UserData');
	nPages    = length(hAxes);
	
	if Cpage~=1
		set(get(hAxes(Cpage),'Children'),'Visible','off'), end
	for p = 1:nPages
		set(get(hAxes(p),'Children'),'Visible','on')
		eval(PrintCmd)
		set(get(hAxes(p),'Children'),'Visible','off')
	end
	set(get(hAxes(Cpage),'Children'),'Visible','on')
end
return


elseif strcmp(Action,'CreateBar')
%=======================================================================
% spm_figure('CreateBar',F)

if (nargin<2)
	if any(get(0,'Children'))
		F = gcf;
	else
		error('No figures available to create tool bar in!')
	end
else
	F = P2;
end

F = spm_figure('FindWin',F);
if isempty(F), error('Figure not found'), end

%-Get position and size parameters
%-----------------------------------------------------------------------
figure(F);
set(F,'Units','Pixels');
P     = get(F,'Position'); P  = P(3:4);		% Figure dimensions {pixels}
S_Gra = P./[600, 865];				% x & y scaling coefs

nBut  = 10;
nGap  = 3;
sx    = floor(P(1)./(nBut+(nGap+2)*0.4));	% uicontrol object width
dx    = floor(2*sx/5);				% inter-uicontrol gap
sy    = floor(20*S_Gra(1));			% uicontrol object height
x     = dx;					% initial x position
y     = P(2) - sy;				% uicontrol y position

%-Delete any existing uicontrol objects
%-----------------------------------------------------------------------
h = findobj(F,'Tag','NoDelete');
delete(h);

%-Create uicontrol objects
%-----------------------------------------------------------------------
uicontrol(F,'Style', 'Frame',...
	'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
	'Tag','NoDelete');

uicontrol(F,'String','Print' ,'Position',[x y sx sy],...
	'CallBack','spm_figure(''Print'',gcf)',...
	'Interruptible','No',...
	'Tag','NoDelete','ForegroundColor','b'); x = x+sx+dx;

h = uicontrol(F,'String','Clear' ,'Position',[x y sx sy],...
	'CallBack','spm_figure(''Clear'',gcf)',...
	'Interruptible','No',...
        'Tag','NoDelete','ForegroundColor','b'); x = x+sx+dx;
if strcmp(get(F,'Tag'),'Graphics')
	set(h,'CallBack',...
	'spm_figure(''Clear'',gcf), spm_figure(''Clear'',''Interactive'')')
end

uicontrol(F,'Style','PopUp','String','gray|hot|split',...
	'Position',[x,y,2*sx,sy],...
	'CallBack','spm_figure(''ColorMap'',get(gco,''Value''))',...
	'Tag','NoDelete',...
	'UserData','Pop'); x = x + 2*sx;

uicontrol(F,'Style','PopUp','String','invert|brighten|darken',...
	'Position',[x,y,2*sx,sy],...
	'CallBack','spm_figure(''ColorMap'',get(gco,''Value''))',...
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

return


elseif strcmp(Action,'ColorMap')
%=======================================================================
% spm_figure('Colormap',ColAction)
if nargin<2, ColAction='gray'; else, ColAction=P2; end
if ~isstr(ColAction)
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
elseif strcmp(ColAction, 'brighten')
	colormap(brighten(colormap, 0.2))
elseif strcmp(ColAction, 'darken')
	colormap(brighten(colormap, -0.2))
else
	error('Illegal ColAction specification')
end

return

elseif strcmp(Action,'GraphicsCut')
%=======================================================================
% Delete next object clicked, provided it's deletable and in this figure
% spm_figure('GraphicsCut')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
waitforbuttonpress;

h = gco(F);
NoDel=[F; findobj(F,'Tag','NoDelete')];
if ( ~any(h==NoDel) & (gcf==F) ), delete(h); end

set(hBut,'ForegroundColor','k')

return


elseif strcmp(Action,'GraphicsMoveStart')
%=======================================================================
% Move the next object clicked, provided it's movable and in this figure
% spm_figure('GraphicsMove')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
waitforbuttonpress;

hPress = gco(F);
NoDel=[F; findobj(F,'Tag','NoDelete')];
if (~any(hPress==NoDel) & gcf==F)
	set(F,'Units','Pixels')
	OPt  = get(F,'CurrentPoint');
	tmp  = get(hPress,'Type');
	if ( strcmp(tmp,'patch') | strcmp(tmp,'line') | strcmp(tmp,'image') )
		hMove = get(hPress,'Parent');
	else
		hMove = hPress;
	end
	set(hMove,'Units','Pixels');
	OPos = get(hMove,'Position');
	set(hPress,'UserData',[hMove,OPt,OPos])
	%set(F,'UserData',[h,OPt,OPos])
	set(F,'WindowButtonUpFcn','spm_figure(''GraphicsMoveEnd'')')
	set(F,'WindowButtonMotionFcn','spm_figure(''GraphicsMoveMotion'')')
end

set(hBut,'ForegroundColor','k')

return


elseif strcmp(Action,'GraphicsMoveMotion')
%=======================================================================
F     = gcf;
hPress = gco;
tmp   = get(hPress,'UserData');
hMove = tmp(1); OPt = tmp(2:3); OPos = tmp(4:length(tmp));
CPt   = get(F,'CurrentPoint');
Disp  = CPt - OPt;
set(hMove,'Units','Pixels')
CPos  = OPos; CPos(1:2) = CPos(1:2) + Disp;
set(hMove,'Position',CPos)
set(hMove,'Units','Normalized')


elseif strcmp(Action,'GraphicsMoveEnd')
%=======================================================================

set(gcf,'UserData',[])
set(gcf,'WindowButtonMotionFcn',' ')
set(gcf,'WindowButtonUpFcn',' ')

return



elseif strcmp(Action,'GraphicsSize')
%=======================================================================
% Change size of next object clicked, provided it's editable and in figure
% spm_figure('GraphicsSize')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
waitforbuttonpress;

h = gco(F);
NoDel=[F; findobj(F,'Tag','NoDelete')];
if ( ~any(h==NoDel) & (gcf==F) )
	u = strcmp(get(F,'SelectionType'),'alt');
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



elseif strcmp(Action,'GraphicsText')
%=======================================================================
% Add text annotation to a figure
% spm_figure('GraphicsText')

F = gcf;
hBut  = gco;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name','Select starting position, then edit text widget...');
set(F,'Pointer','BotL')
waitforbuttonpress;
set(F,'Pointer','Arrow')
set(F,'Name',tmp)

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



else
%=======================================================================
error('Illegal Action string')


%=======================================================================
end
