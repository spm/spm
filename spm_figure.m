function R1=spm_figure(Action,P2,P3,P4)
% Setup and callback function for graphics window
% FORMAT R1=spm_figure(Action,P2,P3,P4)
%	- An embedded callback, multi-function function
%
% FORMAT F = spm_figure
% Defaults to Action 'Create'
%
% FORMAT F = spm_figure(F)
% [ShortCut] Where F is numeric, defaults to spm_figure('Create',F)
%
% FORMAT F = spm_figure('Create',F)
% Create a ToolBar in figure with 'Tag' string or figure number F.
% If F is omitted then the current figure is used, or a new one
% created (spm_figure('CreateWin')).
% F	- (Input)  Figure to use [Optional] - 'Tag' string or figure number.
% F	- (Output) Figure used
%
% FORMAT F = spm_figure('Create')
% Creates a full length WhiteBg figure tagged 'Graphics' for results.
% F	- (Output) Figure created
%
% FORMAT F = spm_figure('FindWin',F)
% Finds window with 'Tag' or figure numnber F - returns empty F if not found
% F	- (Input)  Figure to use [Optional] - 'Tag' string or figure number.
%	- Defaults to 'Graphics'
% F	- (Output) Figure number (if found) or empty (if not).
%
% FORMAT spm_figure('Clear',F)
% Clears figure, redraws ToolBar, clears linked windows specified
% in 'UserData' of figure F.
% F	- 'Tag' string or figure number of figure to clear
%
% FORMAT spm_figure('CreateBar',F)
% Creates toolbar in figure F
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
%
% spm_figure creates user interface objects in the 'results window'
% that faciliate interactive editing of graphic output prior to
% printing (e.g. selection of color maps, deleting, moving and
% editing graphics objects or adding text)
%
% Objects with the attribute Tag = 'NoDelete' are exempt from deletion
% when spm_clf is used.
%
% For SPM usage, the figure should be 'Tag@ged as 'Graphics', in which
% case "Clear" will clear the 'Interactive' window as well.
%
% see also: spm_print
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

if (nargin==0), Action = 'Create'; end
if ~isstr(Action), P2=Action; Action='Create'; end

if strcmp(Action,'Create')
%=======================================================================
% F = spm_figure('Create',F)
%-Condition arguments
if (nargin<2), F=[]; else, F=P2; end

%-Find the window specified. If not specified create new
if isempty(F)
	F = spm_figure('CreateWin');
else
	%-Check out the specified figure
	F = spm_figure('FindWin',F);
	if isempty(F), error('Invalid figure handle'), end
end

spm_figure('CreateBar',F)

R1 = F;

return



elseif strcmp(Action,'CreateWin')
%=======================================================================
S0     = get(0,'ScreenSize');
Fsca   = [S0(3)/1152 S0(4)/900]; Fsca = [Fsca,Fsca];
S_Gra  = [515 008 600 865].*Fsca;
F      = figure('Name','Results',...
	'NumberTitle','off',...
	'Position',S_Gra,...
	'Resize','off',...
	'Visible','on',...
	'PaperPosition',[.75 1.5 7 9.5]);
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
%-Returns -1 if window cannot be found - deletes multiple tagged figs.

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
%-Clear Graphics window, reset buttons, clear linked windows too.

if (nargin<2)
	error('Insufficient arguments: Specify figure')
else
	F = P2;
end

F=spm_figure('FindWin',F);
if isempty(F), error('Invalid figure handle'), end

%-Clear Graphics figure, redraw ToolBar
figure(F), clf, set(F,'Pointer','Arrow')
spm_figure('CreateBar',F);

%-If in SPM usage - clear any 'Interactive' tagged windows.
if strcmp(get(F,'Tag'),'Graphics')
	OF=findobj(get(0,'Children'),'Flat','Tag','Interactive');
	for of = OF(:)'
		figure(of), clf, set(of,'Pointer','Arrow')
	end
end

return



elseif strcmp(Action,'CreateBar')
%=======================================================================
% spm_figure('CreateBar',F)

if (nargin<2)
	error('Insufficient arguments: Figure # required')
else
	F = P2;
end

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

%-Create uicontrol objects
%-----------------------------------------------------------------------
uicontrol(F,'Style', 'Frame',...
	'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
	'Tag','NoDelete');

uicontrol(F,'String','Print' ,'Position',[x y sx sy],...
	'CallBack','spm_print(gcf)',...
	'Interruptible','No',...
	'Tag','NoDelete','ForegroundColor','b'); x = x+sx+dx;

uicontrol(F,'String','Clear' ,'Position',[x y sx sy],...
	'CallBack','spm_figure(''Clear'',gcf)',...
	'Interruptible','No',...
        'Tag','NoDelete','ForegroundColor','b'); x = x+sx+dx;

uicontrol(F,'String','gray'  ,'Position',[x y sx sy],...
	'CallBack','colormap(gray(64))',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','hot'   ,'Position',[x y sx sy],...
	'CallBack','colormap(hot(64))' ,...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','split' ,'Position',[x y sx sy],...
	'CallBack','load Split; colormap(split)',...
	'Tag','NoDelete'); x = x+sx;
uicontrol(F,'String','invert','Position',[x y sx sy],...
	'CallBack','colormap(flipud(colormap))',...
	'Tag','NoDelete'); x = x+sx+dx;

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
		% spm_position(h,u,[0 0]);
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
