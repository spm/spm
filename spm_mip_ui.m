function varargout=spm_mip_ui(varargin)
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
% FORMAT spm_mip_ui('HideGreens',cbf)
% Hides green secondary marker points
% cbf     - [Optional] If true then resets WindowButtonUpFcn
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

%-Axis offsets for 3d MIPs: 
%-----------------------------------------------------------------------
%-MIP pane dimensions and Talairach origin offsets
%-See spm_project.c for derivation
DXYZ = [182 218 182];
CXYZ = [091 127 073];
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
Po(1)  =                  CXYZ(2) -4;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1)   ;
Po(3)  = DXYZ(3)         -CXYZ(3)   ;
Po(4)  = DXYZ(2)         +CXYZ(1) -4;


%-Condition arguments
%=======================================================================
if nargin==0
	error('Insufficient arguments')
elseif ~isstr(varargin{1})
	spm_mip_ui('Display',varargin{1:end}), return
end


switch lower(varargin{1}), case 'display'
%=======================================================================
% spm_mip_ui('Display',t,XYZ,V,Descrip)
if nargin < 4, error('Insufficient arguments'), end
t       = varargin{2};
XYZ     = varargin{3};
V       = varargin{4};
if nargin < 5, Descrip='SPM Maximum Intensity Projections';
	else, Descrip = varargin{5}; end


%-Display MIP
%-----------------------------------------------------------------------
figure(Fgraph)
spm_clf(Fgraph)
set(Fgraph,'Units','normalized')
hMIPax = axes('Position',[0.24 0.54 0.62 0.42]);
spm_mip(t,XYZ,V);

%-Canonicalise axis positioning & save
%-----------------------------------------------------------------------
%-We're assumming the MIP fills the MIP axes bounding box, which is only
% true if stretch-to-fit is on (which it isn't with `axis image`), or that
% the axis bounding box itself has the correct aspect ratio.
%-So, shrink axes to the correct AspectRatio.
% Get aspect ratio from MIP graphic in case 2D
AR       = diff(get(hMIPax,'Xlim'))/diff(get(hMIPax,'Ylim'));
set(hMIPax,'Units','pixels')
MIPaxPos = get(hMIPax,'Position');
if (MIPaxPos(3)/MIPaxPos(4) > AR), MIPaxPos(3) = MIPaxPos(4)*AR;
	else, MIPaxPos(4) = MIPaxPos(3)./AR; end
set(hMIPax,'Position',MIPaxPos,'Tag','MIPaxes','UserData',MIPaxPos(1:2))

%-Create point markers
%-----------------------------------------------------------------------
xyz = spm_mip_ui('RoundCoords',[0,0,0],V);
if V(3) == 1
	%-2 dimensional data
	%---------------------------------------------------------------
	X1   = text(xyz(1),xyz(2),'<','Color','r','Fontsize',16,...
			'Tag','X1',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X1g  = text(xyz(1),xyz(2),'<','Color','g','Fontsize',12,...
			'Tag','X1g','Visible','off');
else
	%-Create point markers
	%--------------------------------------------------------------
	X1   = text(Po(1)+xyz(2),Po(2)+xyz(1),'<',...
			'Color','r','Fontsize',16,...
			'Tag','X1',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X2   = text(Po(1)+xyz(2),Po(3)-xyz(3),'<',...
			'Color','r','Fontsize',16,...
			'Tag','X2',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	X3   = text(Po(4)+xyz(1),Po(3)-xyz(3),'<',....
			'Color','r','Fontsize',16,...
			'Tag','X3',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');

	X1g  = text(Po(1),Po(2),'<','Color','g','Fontsize',12,...
			'Tag','X1g','Visible','off');
	X2g  = text(Po(1),Po(3),'<','Color','g','Fontsize',12,...
			'Tag','X2g','Visible','off');
	X3g  = text(Po(4),Po(3),'<','Color','g','Fontsize',12,...
			'Tag','X3g','Visible','off');
end

%-Print coordinates
%-----------------------------------------------------------------------
hCoordAxes = axes('Position',[0.08 0.80 0.10 0.03],...
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
	'Position',[025 745 090 016].*WS,...
	'FontName','Times','FontWeight','Bold',...
	'ForegroundColor','w',...
	'HorizontalAlignment','Center')
uicontrol(Fgraph,'Style','Text','FontName','Times',...
	'Position',[30 720 16 16].*WS,'String','x');
uicontrol(Fgraph,'Style','Text','FontName','Times',...
	'Position',[30 700 16 16].*WS,'String','y');
uicontrol(Fgraph,'Style','Text','FontName','Times',...
	'Position',[30 680 16 16].*WS,'String','z');

%-Editable co-ordinate windows
%-----------------------------------------------------------------------
hX   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(1)),...
	'FontSize',12,...
	'Tag','hX',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 718 055 020].*WS);
hY   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(2)),...
	'FontSize',12,...
	'Tag','hY',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 698 055 020].*WS);
hZ   = uicontrol(Fgraph,'Style','Edit',...
	'String',sprintf('%.2f',xyz(3)),...
	'FontSize',12,...
	'Tag','hZ',...
	'Callback','spm_mip_ui(''SetCoords'');',...
	'Position',[050 678 055 020].*WS);


case 'getcoords'
%=======================================================================
% xyz = spm_mip_ui('GetCoords')

hX = findobj(Fgraph,'Tag','hX');
hY = findobj(Fgraph,'Tag','hY');
hZ = findobj(Fgraph,'Tag','hZ');

if length([hX,hY,hZ])~=3, error('Co-ord widgets not found'), end

xyz = [	eval(get(hX,'String'));...
	eval(get(hY,'String'));...
	eval(get(hZ,'String'))	];

varargout = {xyz};


case 'roundcoords'
%=======================================================================
% xyz = spm_mip_ui('RoundCoords',xyz,V)
if nargin < 3, V = []; else, V = varargin{3}; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = varargin{2}; end

%-Round xyz to coordinates of actual voxel centre
%-Do rounding in voxel coordinates & ensure within image size
%-----------------------------------------------------------------------
M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
xyz = [xyz(:); 1];
rcp = round(inv(M)*xyz);
rcp = max([min([rcp';[V(1:3)',1]]);[1,1,1,1]])';
xyz = M*rcp;

varargout = {xyz(1:3)};


case 'setcoords'
%=======================================================================
% xyz = spm_mip_ui('SetCoords',xyz,V)
if nargin < 3, V = []; else, V = varargin{3}; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = varargin{2}; end

xyz = spm_mip_ui('RoundCoords',xyz,V);

spm_mip_ui('PrntCords',xyz)
spm_mip_ui('SetCordCntls',xyz)
spm_mip_ui('PosnMarkerPoints',xyz,V)

varargout = {xyz(1:3)};


case 'prntcords'
%=======================================================================
% spm_mip_ui('PrntCords',xyz)
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = varargin{2}; end

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


case 'setcordcntls'
%=======================================================================
% spm_mip_ui('SetCordCntls',xyz)
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = varargin{2}; end

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


case 'posnmarkerpoints'
%=======================================================================
% spm_mip_ui('PosnMarkerPoints',xyz,V)
if nargin < 3, V = []; else, V = varargin{3}; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = spm_mip_ui('GetCoords'); else, xyz = varargin{2}; end

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
	set(X1,'Position',[ Po(1) + xyz(2), Po(2) + xyz(1), 0])
	set(X2,'Position',[ Po(1) + xyz(2), Po(3) - xyz(3), 0])
	set(X3,'Position',[ Po(4) + xyz(1), Po(3) - xyz(3), 0])
end


case 'showgreens'
%=======================================================================
% spm_mip_ui('ShowGreens',xyz,V)
if nargin < 3, V = []; else, V = varargin{3}; end
if isempty(V), V = get(findobj('Tag','hV'),'UserData'); end
if nargin < 2, xyz = get(gco,'UserData'); else, xyz = varargin{2}; end

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
	set(X1g,'Position',[ Po(1) + xyz(2), Po(2) + xyz(1), 0],'Visible','on')
	set(X2g,'Position',[ Po(1) + xyz(2), Po(3) - xyz(3), 0],'Visible','on')
	set(X3g,'Position',[ Po(4) + xyz(1), Po(3) - xyz(3), 0],'Visible','on')
end

%-If called as a callback (xyz from gco), then set WindowButtonUpFcn
%-----------------------------------------------------------------------
if nargin < 2
	set(gcbf,'WindowButtonUpFcn','spm_mip_ui(''HideGreens'',1)')
end


case 'hidegreens'
%=======================================================================
% spm_mip_ui('HideGreens',cbf)
if nargin<2, cbf=0; else, cbf=1; end

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

%-Reset WindowButtonUpFcn if this was a callback function
%-----------------------------------------------------------------------
if cbf, set(gcbf,'WindowButtonUpFcn',' '), end


case 'movestart'
%=======================================================================
% spm_mip_ui('MoveStart')
[cO,cF] = gcbo;

%-Store useful quantities in UserData of gcbo, the object to be dragged
%-----------------------------------------------------------------------
MS.xyz      = spm_mip_ui('GetCoords');
MS.MIPaxPos = get(findobj(cF,'Tag','MIPaxes'),'UserData');
MS.V        = get(findobj(cF,'Tag','hV'),'UserData');
MS.X1       = findobj(cF,'Tag','X1');
if MS.V(3) == 1
	MS.X2 = 0;
	MS.X3 = 0;
else
	MS.X2 = findobj(cF,'Tag','X2');
	MS.X3 = findobj(cF,'Tag','X3');
end
MS.hX      = findobj(cF,'Tag','hX');
MS.hY      = findobj(cF,'Tag','hY');
MS.hZ      = findobj(cF,'Tag','hZ');
set(cO,'UserData',MS)

%-Initiate dragging
%-----------------------------------------------------------------------
if strcmp(get(cF,'SelectionType'),'normal')
	%-Set Figure callbacks for drag'n'drop (DragType 1)
	%---------------------------------------------------------------
	set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',1)',...
		'Interruptible','off')
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
		'Interruptible','off')
	set(cF,'Pointer','Fleur')
elseif strcmp(get(cF,'SelectionType'),'extend')
	%-Set Figure callbacks for drag'n'drop with co-ord updating (DragType 2)
	%---------------------------------------------------------------
	set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',2)',...
		'Interruptible','off')
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'',2)',...
		'Interruptible','off')
	set(cF,'Pointer','Fleur')
elseif strcmp(get(cF,'SelectionType'),'alt')
	%-Set Figure callbacks for drop but no drag (DragType 0)
	%---------------------------------------------------------------
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
		'Interruptible','off')
	set(cF,'Pointer','CrossHair')
end


case 'move'
%=======================================================================
% spm_mip_ui('Move',DragType)
if nargin<2, DragType = 2; else, DragType = varargin{2}; end
cF = gcbf;
cO = gco;

%-Get useful data from UserData of gcbo, the object to be dragged
%-----------------------------------------------------------------------
MS = get(cO,'UserData');
b2d = MS.V(3) == 1;

%-Work out where we are moving to - Use HandleGraphics to give positon
%-----------------------------------------------------------------------
set(cF,'Units','pixels')
d = get(cF,'CurrentPoint') - MS.MIPaxPos;
set(cO,'Units','pixels')
set(cO,'Position',d)
set(cO,'Units','data')
d = get(cO,'Position');


%-Work out xyz, depending on which view is being manipulated
%-----------------------------------------------------------------------
sMarker = get(cO,'Tag');
if strcmp(sMarker,'X1')
	if b2d
		xyz = [d(1); d(2); MS.xyz(3)];
	else
		xyz = [d(2) - Po(2); d(1) - Po(1); MS.xyz(3)];
	end
elseif strcmp(sMarker,'X2')
	xyz = [MS.xyz(1); d(1) - Po(1); Po(3) - d(2)];
elseif strcmp(sMarker,'X3')
	xyz = [d(1) - Po(4); MS.xyz(2); Po(3) - d(2)];
else
	error('Can''t work out which marker point')
end

%-Round coordinates to nearest voxel center
%-----------------------------------------------------------------------
% xyz = spm_mip_ui('RoundCoords',xyz,V);
M = [ [diag(MS.V(4:6)), -(MS.V(7:9).*MS.V(4:6))]; [zeros(1,3) ,1]];
xyz = M*max([min([round(inv(M)*[xyz(:); 1])';[MS.V(1:3)',1]]);[1,1,1,1]])';
xyz = xyz(1:3);

%-Move marker points
%-----------------------------------------------------------------------
% spm_mip_ui('PosnMarkerPoints',xyz,MS.V)
if b2d
	set(MS.X1,'Units','Data')
	set(MS.X1,'Position',[xyz(1), xyz(2), 1],'Visible','on')
else
	set([MS.X1,MS.X2,MS.X3],'Units','Data')
	set(MS.X1,'Position',[ Po(1) + xyz(2), Po(2) + xyz(1), 0])
	set(MS.X2,'Position',[ Po(1) + xyz(2), Po(3) - xyz(3), 0])
	set(MS.X3,'Position',[ Po(4) + xyz(1), Po(3) - xyz(3), 0])
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
	set(MS.hX,'String',sprintf('%0.2f',xyz(1)))
	set(MS.hY,'String',sprintf('%0.2f',xyz(2)))
	set(MS.hZ,'String',sprintf('%0.2f',xyz(3)))
else
	error('Illegal DragType')
end


case 'moveend'
%=======================================================================
% spm_mip_ui('MoveEnd',EndType)
if nargin<2, EndType = 0; else, EndType = varargin{2}; end

%-Reset WindowButton functions & pointer
%-----------------------------------------------------------------------
set(gcf,'WindowButtonMotionFcn',' ')
set(gcf,'WindowButtonUpFcn',' ')
set(gcf,'Pointer','arrow')

%-Print coordinates after drag'n'drop
%-----------------------------------------------------------------------
if EndType, spm_mip_ui('PrntCords'), end


case 'showhiddencoords'
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


otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end
