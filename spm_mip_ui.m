function R1=spm_mip_ui(Action,P2,P3,P4,P5)
% UI for displaying MIPs with interactive controls
% FORMAT spm_mip_ui(t,XYZ,V,Descrip)
% t       - SPM point list for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% V       - {9 x 1} vector of image and voxel sizes, and origin
%               [DIM,VOX,ORIGIN]'
% Descrip - [Optional] string matrix of text describing the MIP
%           The first line is printed in bold.
%
% FORMAT xyz = spm_mip_ui('GetCoords')
% xyz     - (Input) {3 x 1} vector of Talairach coordinates for cursor
% xyz     - (Output) Talairach coordinates of cursor, at centre of voxel
%           nearest specified coordinates
%
% FORMAT xyz = spm_mip_ui('SetCoords',xyz,V)
% xyz     - {3 x 1} vector of Talairach coordinates
% V       - [Optional] {9 x 1} vector of image and voxel sizes, and origin
%_______________________________________________________________________
%
% spm_mip_ui displays a maximum intensity projection (using spm_mip)
% with draggable cursors and editable coordinate windows.
% 
% The cursor can be dragged, or the coordinates entered in the editable
% coordinate windows. Dragging with the primary "select" mouse button
% moves the cursor across voxel centers. Dragging with the "extend"
% mouse button simultaneously updates coordinates and cursor position.
% If either are too slow then dragging with the "alt" mouse button only
% updates the cursor positions and coordinates when the mouse button is
% released.
% 
% The current cursor position (constrained to lie on a voxel) can be
% obtained by xyz=spm_mip_ui('GetCoords'), and set with
% xyz=spm_mip_ui('SetCoords',xyz). The latter rounds xyz to the nearest
% voxel center, returning the result.
% 
% spm_mip_ui handles all the callbacks required for moving the cursors.
% Programers help is below.
% 
% The MIP is displayed in the 'Graphics' 'Tag'ged figure window. If
% this does not exist then the current figure is 'Tag'ged 'Graphics',
% or an SPM Graphics
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
% FORMAT spm_mip_ui(t,XYZ,V,Descrip)
% [ShortCut] Defaults to spm_mip_ui('Display',t,XYZ,V,Descrip)
%
% FORMAT spm_mip_ui('Display',t,XYZ,V,Descrip)
% Displays the MIP and sets up cursors and control objects
% t       - SPM point list for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% V       - {9 x 1} vector of image and voxel sizes, and origin
%               [DIM,VOX,ORIGIN]'
% Descrip - [Optional] string matrix of text describing the MIP
%           The first line is printed in bold.
%
% FORMAT xyz = spm_mip_ui('GetCoords')
% Returns coordinates of current cursor position
% xyz     - (Input) {3 x 1} vector of Talairach coordinates for cursor
% xyz     - (Output) Talairach coordinates of cursor, at centre of voxel
%           nearest specified coordinates
%
% FORMAT xyz = spm_mip_ui('RoundCoords',xyz,V)
% Rounds coordinates to centre of nearest voxel
% xyz     - (Input) {3 x 1} vector of Talairach coordinates for cursor
% xyz     - (Output) Talairach coordinates of cursor, at centre of voxel
%           nearest specified coordinates
% V       - [Optional] {9 x 1} vector of image and voxel sizes, and origin
%
%
% FORMAT xyz = spm_mip_ui('SetCoords',xyz,V)
% Sets cursor position
% xyz     - {3 x 1} vector of Talairach coordinates
% V       - [Optional] {9 x 1} vector of image and voxel sizes, and origin
%
% FORMAT spm_mip_ui('PrntCords',xyz)
% Utility routine: Updates printed coordinates
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
%
% FORMAT spm_mip_ui('SetCordCntls',xyz)
% Utility routine: Sets values of coordinate controls
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
%
% FORMAT spm_mip_ui('PosnMarkerPoints',xyz,V)
% Utility routine: Positions cursor markers
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
% V       - [Optional] {9 x 1} vector of image and voxel sizes, and origin
%
% FORMAT spm_mip_ui('ShowGreens',xyz,V)
% Shows green secondary cursors at location xyz
% xyz     - [Optional] {3 x 1} vector of Talairach coordinates for cursor
%            Defaults to UserData of gco
% V       - [Optional] {9 x 1} vector of image and voxel sizes, and origin
%
% FORMAT spm_mip_ui('HideGreens',Win)
% Hides green secondary marker points
% Win     - [Optional] If true then resets WindowButtonUpFcn
%           Used in CallBacks
%
% FORMAT spm_mip_ui('MoveStart')
% Utility routine: CallBack for starting cursor dragging
% This is the ButtonDownFcn for the cursor markers
%
% FORMAT spm_mip_ui('Move',DragType)
% Utility routine: Initiate cursor move
% DragType - 0 = Drop (no dragging)
%            1 = Drag'n'drop
%            2 = Drag'n'drop with dynamic coordinate updating
%
% FORMAT spm_mip_ui('MoveEnd',EndType)
% Utility routine: End cursor move
% EndType - Flag. If set then MoveEnd updates printed coordinates
%
% FORMAT spm_mip_ui('ShowHiddenCoords')
% Utility routine: Moves printed coordinate axes into/out of view
% ( Printed coordinate axes are hidden behind the coordinate controls. This  )
% ( is required since the uicontrols are invisible on printouts.             )
% This function is the callbace to the first "Descrip"tion line.
% The "extend" mouse button toggles between show/hide. The "alt" button
% shows the text axes while the button is depressed.
%
%_______________________________________________________________________


%-Global Parameters
%=======================================================================
WS     = spm('GetWinScale');
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph)
	if any(get(0,'Children')), F = gcf; set(F,'Tag','Graphics')
		else, Fgraph=spm_figure('Create','Graphics'); end
end

%-Axis offsets for 3d MIPs: See spm_project.c for derivation
%-----------------------------------------------------------------------
Po1 = 127-4;
Po2 = 182+182-91;
Po3 = 182-73;
Po4 = 218+91-4;


%-Format arguments
%=======================================================================
if nargin == 0, Action='Welcome'; end
if ~isstr(Action) & nargin>=3
	if nargin==4, P5 = P4; end
	P4 = P3;
	P3 = P2;
	P2 = Action;
	Action = 'Display';
	nargin = nargin + 1;
end


if strcmp(lower(Action),lower('Display'))
%=======================================================================
% spm_mip_ui('Display',t,XYZ,V,Descrip)
if nargin < 4, error('Insufficient arguments'), end
t       = P2;
XYZ     = P3;
V       = P4;
if nargin < 5, Descrip='SPM Maximum Intensity Projections';
	else, Descrip = P5; end

%-Display MIP
%-----------------------------------------------------------------------
figure(Fgraph)
spm_clf(Fgraph)
set(Fgraph,'Units','normalized'); 
hMIPaxes = axes('Position',[0.24 0.54 0.62 0.42],'Tag','hMIPaxes');
spm_mip(t,XYZ,V); axis normal


%-Create point markers
%-----------------------------------------------------------------------
xyz = spm_mip_ui('RoundCoords',[0,0,0],V);
if V(3) == 1
	%-2 dimensional data
	%---------------------------------------------------------------
	set(hMIPaxes,'Position',[0.24 0.54 0.52 0.36])
	X1   = text(xyz(1),xyz(2),'<','Color','r','Fontsize',16,...
			'Tag','X1',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X1g  = text(xyz(1),xyz(2),'<','Color','g','Fontsize',12,...
			'Tag','X1g','Visible','off');
else
	%-Create point markers
	%--------------------------------------------------------------
	X1   = text(Po1+xyz(2),Po2+xyz(1),'<','Color','r','Fontsize',16,...
			'Tag','X1',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X2   = text(Po1+xyz(2),Po3-xyz(3),'<','Color','r','Fontsize',16,...
			'Tag','X2',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X3   = text(Po4+xyz(1),Po3-xyz(3),'<','Color','r','Fontsize',16,...
			'Tag','X3',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');

	X1g  = text(Po1,Po2,'<','Color','g','Fontsize',12,...
			'Tag','X1g','Visible','off');
	X2g  = text(Po1,Po3,'<','Color','g','Fontsize',12,...
			'Tag','X2g','Visible','off');
	X3g  = text(Po4,Po3,'<','Color','g','Fontsize',12,...
			'Tag','X3g','Visible','off');

end

%-Reset axis attributes and get hMIPaxes position in pixels
%-----------------------------------------------------------------------
set(hMIPaxes,'Units','Pixels')
set(Fgraph,'Units','Pixels')
AXES = get(hMIPaxes,'Position');
AXES = AXES(1:2);
set(hMIPaxes,'UserData',AXES,...
	'Units','normalized')

%-Print coordinates
%-----------------------------------------------------------------------
% tmp = [0.08 0.90 0.10 0.03];
tmp = [0.08 0.80 0.10 0.03];
hCoordAxes = axes('Position',tmp,...
	'Visible','off',...
	'Units','Normalized',...
	'Tag','hV',...
	'UserData',V);

text(-0.1,-0.5,'SPM - MIP',...
	'FontSize',8,...
	'FontName','Times','FontWeight','Bold','FontAngle','Oblique',...
	'Rotation',90)

hXstr = text(0,1.0,sprintf('x = %.2f',xyz(1)),'FontSize',10,...
		'Tag','hXstr');
hYstr = text(0,0.5,sprintf('y = %.2f',xyz(2)),'FontSize',10,...
		'Tag','hYstr');
hZstr = text(0,0.0,sprintf('z = %.2f',xyz(3)),'FontSize',10,...
		'Tag','hZstr');

%-Print parameter text
%-----------------------------------------------------------------------
hTextAxes = axes('Position',[0.625 0.48 0.3 0.15],'Visible','off',...
	'Units','Normalized','Tag','hTextAxes');
h = text(0.0,1.0,deblank(Descrip(1,:)),'FontSize',8,'FontWeight','Bold');
y = 0.8; dy = 0.1;
for i = 2:size(Descrip,1)
	text(0,y,deblank(Descrip(i,:)),'FontSize',8)
	y = y - dy;
end
set(h,'UserData',hCoordAxes,...
	'ButtonDownFcn','spm_mip_ui(''ShowHiddenCoords'')');

%-Create control objects
%=======================================================================

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fgraph,'Style','Frame','Position',[20 670 100 096].*WS);
uicontrol(Fgraph,'Style','Text','String','SPM - MIP',...
	'Position',[035 740 070 016],...
	'ForegroundColor','w',...
	'HorizontalAlignment','Center')
uicontrol(Fgraph,'Style','Text','Position',[30 720 16 16].*WS,'String','x');
uicontrol(Fgraph,'Style','Text','Position',[30 700 16 16].*WS,'String','y');
uicontrol(Fgraph,'Style','Text','Position',[30 680 16 16].*WS,'String','z');

%-Editable co-ordinate windows
%-----------------------------------------------------------------------
hX   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(1)),...
	'Tag','hX',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 720 055 016].*WS);
hY   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(2)),...
	'Tag','hY',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 700 055 016].*WS);
hZ   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(3)),...
	'Tag','hZ',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 680 055 016].*WS);
return


elseif strcmp(lower(Action),lower('GetCoords'))
%=======================================================================
% xyz = spm_mip_ui('GetCoords')

hX = findobj(Fgraph,'Tag','hX');
hY = findobj(Fgraph,'Tag','hY');
hZ = findobj(Fgraph,'Tag','hZ');

xyz = [	eval(get(hX,'String'));...
	eval(get(hY,'String'));...
	eval(get(hZ,'String'))	];

R1 = xyz;
return


elseif strcmp(lower(Action),lower('RoundCoords'))
%=======================================================================
% xyz = spm_mip_ui('RoundCoords',xyz,V)
if nargin < 3, V = []; else, V = P3; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = P2; end

%-Round xyz to coordinates of actual voxel centre
%-Do rounding in voxel coordinates & ensure within image size
%-----------------------------------------------------------------------
M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
xyz = [xyz(:); 1];
rcp = round(inv(M)*xyz);
rcp = max([min([rcp';[V(1:3)',1]]);[1,1,1,1]])';
xyz = M*rcp;

R1 = xyz(1:3);
return


elseif strcmp(lower(Action),lower('SetCoords'))
%=======================================================================
% xyz = spm_mip_ui('SetCoords',xyz,V)
if nargin < 3, V = []; else, V = P3; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = P2; end

xyz = spm_mip_ui('RoundCoords',xyz,V);

spm_mip_ui('PrntCords',xyz)
spm_mip_ui('SetCordCntls',xyz)
spm_mip_ui('PosnMarkerPoints',xyz,V)

R1 = xyz(1:3);
return


elseif strcmp(lower(Action),lower('PrntCords'))
%=======================================================================
% spm_mip_ui('PrntCords',xyz)
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = P2; end

%-Get handles of Co-ordinate strings
%-----------------------------------------------------------------------
hXstr = findobj(Fgraph,'Tag','hXstr');
hYstr = findobj(Fgraph,'Tag','hYstr');
hZstr = findobj(Fgraph,'Tag','hZstr');

%-Set Co-ordinate strings
%-----------------------------------------------------------------------
set(hXstr,'String',sprintf('x = %0.2f',xyz(1)))
set(hYstr,'String',sprintf('y = %0.2f',xyz(2)))
set(hZstr,'String',sprintf('z = %0.2f',xyz(3)))

return


elseif strcmp(lower(Action),lower('SetCordCntls'))
%=======================================================================
% spm_mip_ui('SetCordCntls',xyz)
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = P2; end

%-Get handles of Co-ordinate controls
%-----------------------------------------------------------------------
%-Co-ordinate controls
hX    = findobj(Fgraph,'Tag','hX');
hY    = findobj(Fgraph,'Tag','hY');
hZ    = findobj(Fgraph,'Tag','hZ');

%-Set Co-ordinate control strings
%-----------------------------------------------------------------------
set(hX,'String',sprintf('%0.2f',xyz(1)))
set(hY,'String',sprintf('%0.2f',xyz(2)))
set(hZ,'String',sprintf('%0.2f',xyz(3)))

return


elseif strcmp(lower(Action),lower('PosnMarkerPoints'))
%=======================================================================
% spm_mip_ui('PosnMarkerPoints',xyz,V)
if nargin < 3, V = []; else, V = P3; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = P2; end

b2d = V(3) == 1;

%-Get handles of marker points
%-----------------------------------------------------------------------
X1 = findobj(Fgraph,'Tag','X1');
if ~b2d
	X2 = findobj(Fgraph,'Tag','X2');
	X3 = findobj(Fgraph,'Tag','X3');
end

%-Set marker points
%-----------------------------------------------------------------------
if b2d
	set(X1,'Units','Data')
	set(X1,'Position',[xyz(1), xyz(2), 1])
else
	set([X1,X2,X3],'Units','Data')
	set(X1,'Position',[ Po1 + xyz(2), Po2 + xyz(1), 0])
	set(X2,'Position',[ Po1 + xyz(2), Po3 - xyz(3), 0])
	set(X3,'Position',[ Po4 + xyz(1), Po3 - xyz(3), 0])
end
return


elseif strcmp(lower(Action),lower('ShowGreens'))
%=======================================================================
% spm_mip_ui('ShowGreens',xyz,V)
if nargin < 3, V = []; else, V = P3; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = get(gco,'UserData'); else, xyz = P2; end

b2d = V(3) == 1;

%-Get handles of green marker points
%-----------------------------------------------------------------------
X1g = findobj(Fgraph,'Tag','X1g');
if ~b2d
	X2g = findobj(Fgraph,'Tag','X2g');
	X3g = findobj(Fgraph,'Tag','X3g');
end

%-Set marker points
%-----------------------------------------------------------------------
if b2d
	set(X1g,'Units','Data')
	set(X1g,'Position',[xyz(1), xyz(2), 1],'Visible','on')
else
	set([X1g,X2g,X3g],'Units','Data')
	set(X1g,'Position',[ Po1 + xyz(2), Po2 + xyz(1), 0],'Visible','on')
	set(X2g,'Position',[ Po1 + xyz(2), Po3 - xyz(3), 0],'Visible','on')
	set(X3g,'Position',[ Po4 + xyz(1), Po3 - xyz(3), 0],'Visible','on')
end

%-If called as a callback (xyz from gco), then set WindowButtonUpFcn
%-----------------------------------------------------------------------
if nargin<2
	set(gcf,'WindowButtonUpFcn','spm_mip_ui(''HideGreens'',1)')
end

return


elseif strcmp(lower(Action),lower('HideGreens'))
%=======================================================================
% spm_mip_ui('HideGreens',Win)
if nargin<2, Win=0; else, Win=1; end

V = get(findobj('Tag','hV'),'UserData');
b2d = V(3) == 1;

%-Get handles of green marker points and hide them
%-----------------------------------------------------------------------
X1g = findobj(Fgraph,'Tag','X1g');
set(X1g,'Visible','off')
if ~b2d
	X2g = findobj(Fgraph,'Tag','X2g');
	X3g = findobj(Fgraph,'Tag','X3g');
	set([X1g,X2g,X3g],'Visible','off')
end

%-Reset WindowButtonUpFcn if needed
%-----------------------------------------------------------------------
if Win, set(gcf,'WindowButtonUpFcn',' '), end
return


elseif strcmp(lower(Action),lower('MoveStart'))
%=======================================================================
% spm_mip_ui('MoveStart')

%-Store useful quantities in UserData of gco, the object to be dragged
%-----------------------------------------------------------------------
xyz  = spm_mip_ui('GetCoords');
AXIS = get(findobj(gcf,'Tag','hMIPaxes'),'UserData');
V    = get(findobj(gcf,'Tag','hV'),'UserData');
X1 = findobj(gcf,'Tag','X1');
b2d = V(3) == 1;
if b2d
	X2 = 0;
	X3 = 0;
else
	X2 = findobj(gcf,'Tag','X2');
	X3 = findobj(gcf,'Tag','X3');
end
hX    = findobj(gcf,'Tag','hX');
hY    = findobj(gcf,'Tag','hY');
hZ    = findobj(gcf,'Tag','hZ');
set(gco,'UserData',[xyz; AXIS'; V; X1; X2; X3; hX; hY; hZ])

%-Initiate dragging
%-----------------------------------------------------------------------
if strcmp(get(gcf,'SelectionType'),'normal')
	%-Set Figure callbacks for drag'n'drop (DragType 1)
	%---------------------------------------------------------------
	set(gcf,'WindowButtonMotionFcn','spm_mip_ui(''Move'',1)',...
		'Interruptible','no')
	set(gcf,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
		'Interruptible','no')
	set(gcf,'Pointer','Fleur')
elseif strcmp(get(gcf,'SelectionType'),'extend')
	%-Set Figure callbacks for drag'n'drop with co-ord updating
	% (DragType 2)
	%---------------------------------------------------------------
	set(gcf,'WindowButtonMotionFcn','spm_mip_ui(''Move'',2)',...
		'Interruptible','no')
	set(gcf,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'',2)',...
		'Interruptible','no')
	set(gcf,'Pointer','Fleur')
elseif strcmp(get(gcf,'SelectionType'),'alt')
	%-Set Figure callbacks for drop but no drag (DragType 0)
	%---------------------------------------------------------------
	set(gcf,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
		'Interruptible','no')
	set(gcf,'Pointer','CrossHair')
end
return


elseif strcmp(lower(Action),lower('Move'))
%=======================================================================
% spm_mip_ui('Move',DragType)
if nargin<2, DragType = 2; else, DragType = P2; end
set(gcf,'Units','pixels')

%-Get useful data from UserData of gco, the object to be dragged
%-----------------------------------------------------------------------
tmp = get(gco,'UserData');
% xyz  = tmp(1:3);
AXIS   = tmp(4:5)';
V      = tmp(6:14);
X1     = tmp(15);
X2     = tmp(16);
X3     = tmp(17);
hX     = tmp(18);
hY     = tmp(19);
hZ     = tmp(20);

b2d = V(3) == 1;

%-Work out where we are moving to
%-----------------------------------------------------------------------
d = get(gcf,'CurrentPoint') - AXIS;
set(gco,'Units','pixels')
set(gco,'Position',d)
set(gco,'Units','data')
d = get(gco,'Position');


%-Work out xyz, depending on which view is being manipulated
%-----------------------------------------------------------------------
sMarker = get(gco,'Tag');
if strcmp(sMarker,'X1')
	if b2d
		xyz = [d(1); d(2); tmp(3)];
	else
		xyz = [d(2) - Po2; d(1) - Po1; tmp(3)];
	end
elseif strcmp(sMarker,'X2')
	xyz = [tmp(1); d(1) - Po1; Po3 - d(2)];
elseif strcmp(sMarker,'X3')
	xyz = [d(1) - Po4; tmp(2); Po3 - d(2)];
else
	error('Can''t work out which marker point')
end

%-Round coordinates to nearest voxel center
%-----------------------------------------------------------------------
% xyz = spm_mip_ui('RoundCoords',xyz,V);
M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
xyz = M*max([min([round(inv(M)*[xyz(:); 1])';[V(1:3)',1]]);[1,1,1,1]])';
xyz = xyz(1:3);

%-Move marker points
%-----------------------------------------------------------------------
% spm_mip_ui('PosnMarkerPoints',xyz,V)
if b2d
	set(X1,'Units','Data')
	set(X1,'Position',[xyz(1), xyz(2), 1],'Visible','on')
else
	set([X1,X2,X3],'Units','Data')
	set(X1,'Position',[ Po1 + xyz(2), Po2 + xyz(1), 0])
	set(X2,'Position',[ Po1 + xyz(2), Po3 - xyz(3), 0])
	set(X3,'Position',[ Po4 + xyz(1), Po3 - xyz(3), 0])
end

%-Set Co-ordinate strings (if appropriate DragType)
%-----------------------------------------------------------------------
% spm_mip_ui('SetCordCntls',xyz)
% spm_mip_ui('PrntCords',xyz)

if DragType==0
	spm_mip_ui('SetCordCntls',xyz)
	spm_mip_ui('PrntCords',xyz)
	spm_mip_ui('MoveEnd')
elseif DragType==1
	% Nothing! DragType 1 ends with ButtonUp of DragType 0 (drop)
elseif DragType==2
	set(hX,'String',sprintf('%0.2f',xyz(1)))
	set(hY,'String',sprintf('%0.2f',xyz(2)))
	set(hZ,'String',sprintf('%0.2f',xyz(3)))
else
	error('Illegal DragType')
end

return


elseif strcmp(lower(Action),lower('MoveEnd'))
%=======================================================================
% spm_mip_ui('MoveEnd',EndType)
if nargin<2, EndType = 0; else, EndType = P2; end

%-Reset WindowButton functions & pointer
%-----------------------------------------------------------------------
set(gcf,'WindowButtonMotionFcn',' ')
set(gcf,'WindowButtonUpFcn',' ')
set(gcf,'Pointer','arrow')

%-Print coordinates after drag'n'drop
%-----------------------------------------------------------------------
if EndType, spm_mip_ui('PrntCords'), end
return

elseif strcmp(lower(Action),lower('ShowHiddenCoords'))
%=======================================================================
% spm_mip_ui('ShowHiddenCoords')
if strcmp(get(gcf,'SelectionType'),'normal'), return, end
h = get(gco,'UserData');
tmp = get(h,'Position');
if tmp(2)==0.9
	%-Co-ordinate axes are showing
	set(h,'Position',[0.08 0.8 0.1 0.03])
	set(gcf,'WindowButtonUpFcn',' ')
else
	%-Co-ordinate axes aren't showing
	set(h,'Position',[0.08 0.9 0.1 0.03])
	if strcmp(get(gcf,'SelectionType'),'alt')
		set(gcf,'WindowButtonUpFcn',...
			'spm_mip_ui(''ShowHiddenCoords'')')
	end
end
return


else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
