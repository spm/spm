function varargout=spm_mip_ui(varargin)
% GUI for displaying MIPs with interactive pointers
% FORMAT hMIPax = spm_mip_ui(Z,XYZ,M,DIM,F)
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
% F       - Figure (or axes) to work in [Defaults to gcf]
% hMIPax  - handle of MIP axes
%
% FORMAT xyz = spm_mip_ui('GetCoords',h)
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% xyz     - Current Talairach coordinates of cursor
%
% FORMAT [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
% xyz     - (Input) {3 x 1} vector of desired Talairach coordinates
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% hC      - Handle of calling object, if used as a callback. [Default 0]
% xyz     - (Output) {3 x 1} vector of voxel centre nearest desired xyz co-ords
% d       - Euclidian distance from desired co-ords & nearest voxel 
%_______________________________________________________________________
%
% spm_mip_ui displays a maximum intensity projection (using spm_mip)
% with draggable cursors.
%
% See spm_mip.m for details of MIP construction, display, and the brain
% outlines used.
%                           ----------------
%
% The cursor can be dragged to new locations in three ways:
%
% (1) Point & drop: Using the primary "select" mouse button, click on a
%     cursor and drag the crosshair which appears to the desired location.
%     On dropping, the cursors jump to the voxel centre nearest the drop
%     site.
%
% (2) Dynamic drag & drop: Using the middle "extend" mouse button, click on
%     a cursor and drag it about. The cursor follows the mouse, jumping to
%     the voxel centre nearest the pointer. A dynamically updating
%     information line appears above the MIP and gives the current
%     co-ordinates. If the current voxel centre is in the XYZ pointlist,
%     then the corresponding image value is also printed.
%
% (3) Magnetic drag & drop: As with "Dynamic drag & drop", except the cursors
%     jump to the voxel centre in the pointlist nearest to the cursor. Use 
%     the right "alt" mouse button for "magnetic drag & drop".
%
% In addition a ContextMenu is provided, giving the option to jump the
% cursors to the nearest suprathreshold voxel, the nearest local
% maxima, or to the global maxima. (Right click on the MIP to bring up
% the ContextMenu.) A message in the MatLab command window describes the 
% jump.
%
%                           ----------------
%
% The current cursor position (constrained to lie on a voxel) can be
% obtained by xyz=spm_mip_ui('GetCoords',hMIPax), and set with
% xyz=spm_mip_ui('SetCoords',xyz,hMIPax), where hMIPax is the handle of
% the MIP axes, or of the figure containing a single MIP [default gcf].
% The latter rounds xyz to the nearest voxel center, returning the
% result.
% 
% spm_mip_ui handles all the callbacks required for moving the cursors, and
% is "registry" enabled (See spm_XYZreg.m). Programmers help is below in the
% main body of the function.
% 
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
% FORMAT hMIPax = spm_mip_ui(Z,XYZ,V,F)
% [ShortCut] Defaults to hMIPax=spm_mip_ui('Display',Z,XYZ,V,F)
%
% FORMAT hMIPax = spm_mip_ui('Display',Z,XYZ,M,DIM,F)
% Displays the MIP and sets up cursors
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
% F       - Handle of figure (or axes) to work in [Defaults to gcf]
% hMIPax  - handle of MIP axes
%
% FORMAT xyz = spm_mip_ui('GetCoords',h)
% Returns coordinates of current cursor position
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% xyz     - Current Talairach coordinates of cursor
%
% FORMAT [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
% Sets cursor position
% xyz     - (Input) {3 x 1} vector of desired Talairach coordinates
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% hC      - Handle of calling object, if used as a callback. [Default 0]
% xyz     - (Output) {3 x 1} vector of voxel centre nearest desired xyz co-ords
% d       - Euclidian distance from desired co-ords & nearest voxel 
%
% FORMAT spm_mip_ui('PosnMarkerPoints',xyz,h,r)
% Utility routine: Positions cursor markers
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% r       - 'r' to move visible red cursors, 'g' to move invisible green cursors
%
%
% FORMAT [xyz,d] = spm_mip_ui('Jump',h,loc)
% Utility routine (CallBack of UIcontextMenu) to jump cursor
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% loc     - String defining jump: 'dntmv' - don't move
%                                 'nrvox' - nearest suprathreshold voxel
%                                 'nrmax' - nearest local maxima
%                                 'glmax' - global maxima
% xyz     - co-ordinates of voxel centre jumped to {3 x 1} vector
% d       - (square) Euclidian distance jumped 
%
% FORMAT spm_mip_ui('ShowGreens',xyz,h)
% Shows green secondary cursors at location xyz
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
%           [Default UserData of gco]
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
%
% FORMAT spm_mip_ui('HideGreens',h)
% Hides green secondary marker points
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
%
% FORMAT hMIPax = spm_mip_ui('FindMIPax',h)
% Looks for / checks MIP axes 'Tag'ged 'hMIPax'... errors if no valid MIP axes
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% hMIPax  - Handle of valid MIP axis found (errors if multiple found)
%
% FORMAT spm_mip_ui('MoveStart')
% Utility routine: CallBack for starting cursor dragging
% This is the ButtonDownFcn for the cursor markers
%
% FORMAT spm_mip_ui('Move',DragType)
% Utility routine: Initiate cursor move
% DragType - 0 = Point & drop (no dragging)
%            1 = Drag'n'drop with dynamic coordinate/value updating
%            2 = Magnetic drag'n'drop with dynamic coordinate updating
%
% FORMAT spm_mip_ui('MoveEnd')
% Utility routine: End cursor move
%_______________________________________________________________________


%-Condition arguments
%=======================================================================
if nargin==0
	error('Insufficient arguments')
elseif ~isstr(varargin{1})
	varargout={spm_mip_ui('Display',varargin{1:end})}; return
end


%-Axis offsets for 3d MIPs: 
%=======================================================================
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
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;



%=======================================================================
switch lower(varargin{1}), case 'display'
%=======================================================================
% hMIPax = spm_mip_ui('Display',Z,XYZ,M,DIM,F)
if nargin<5
	F      = gcf;
	hMIPax = [];
else
	F = varargin{6};
	if isstr(F), F=spm_figure('FindWin',F); end
	if ~ishandle(F), error('Invalid handle'), end
	switch get(F,'Type'), case 'figure'
		hMIPax = [];
	case 'axes'
		hMIPax = F;
		F      = get(hMIPax,'Parent');
	otherwise
		error('F not a figure or axis handle')
	end
end
if nargin<4, error('Insufficient arguments'), end
DIM     = varargin{5};
M       = varargin{4};
XYZ     = varargin{3};
Z       = varargin{2};

xyz = spm_XYZreg('RoundCoords',[0;0;0],M,DIM);


%-Display MIP
%-----------------------------------------------------------------------
Funits = get(F,'Units');
set(F,'Units','normalized')
if isempty(hMIPax)
	hMIPax = axes('Position',[0.24 0.54 0.62 0.42],'Parent',F);
else
	axes(hMIPax), cla reset
end

%-NB: spm_mip's `image` uses a newplot, & screws stuff without the figure.
figure(F)
spm_mip(Z,XYZ,M,DIM);

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
set(hMIPax,'Position',MIPaxPos)
hMIPim = get(gca,'Children');

%-Print coordinates
%-----------------------------------------------------------------------
hMIPxyz = text(0,max(get(hMIPax,'YLim'))/2,...
	{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))},...
	'Interpreter','TeX','FontName',spm_platform('font','times'),...
	'Color',[1,1,1]*.7,...
	'HorizontalAlignment','Center',...
	'VerticalAlignment','Bottom',...
	'Rotation',90,...
	'Tag','hMIPxyz',...
	'UserData',xyz);

%-Create point markers
%-----------------------------------------------------------------------
if DIM(3) == 1
	%-2 dimensional data
	%---------------------------------------------------------------
	hXr    = text(xyz(1),xyz(2),'<','Color','r','Fontsize',16,...
			'Tag','hX1r',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	hXg    = text(xyz(1),xyz(2),'<','Color','g','Fontsize',12,...
			'Tag','hX1g','Visible','off');
else
	%-Create point markers
	%--------------------------------------------------------------
	hX1r   = text(Po(1)+xyz(2),Po(2)+xyz(1),'<',...
			'Color','r','Fontsize',16,...
			'Tag','hX1r',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	hX2r   = text(Po(1)+xyz(2),Po(3)-xyz(3),'<',...
			'Color','r','Fontsize',16,...
			'Tag','hX2r',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	hX3r   = text(Po(4)+xyz(1),Po(3)-xyz(3),'<',....
			'Color','r','Fontsize',16,...
			'Tag','hX3r',...
			'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
	hXr = [hX1r,hX2r,hX3r];

	hX1g  = text(Po(1),Po(2),'<','Color','g','Fontsize',12,...
			'Tag','hX1g','Visible','off');
	hX2g  = text(Po(1),Po(3),'<','Color','g','Fontsize',12,...
			'Tag','hX2g','Visible','off');
	hX3g  = text(Po(4),Po(3),'<','Color','g','Fontsize',12,...
			'Tag','hX3g','Visible','off');
	hXg = [hX1g,hX2g,hX3g];
end

%-Create UIContextMenu for marker jumping
%-----------------------------------------------------------------------
h = uicontextmenu('Tag','MIPconmen','UserData',hMIPax);
uimenu(h,'Label','MIP')
if isempty(XYZ), str='off'; else, str='on'; end
uimenu(h,'Separator','on','Label','goto nearest suprathreshold voxel',...
	'CallBack',['spm_mip_ui(''Jump'',',...
		'get(get(gcbo,''Parent''),''UserData''),''nrvox'');'],...
	'Interruptible','off','BusyAction','Cancel','Enable',str);
uimenu(h,'Separator','off','Label','goto nearest local maxima',...
	'CallBack',['spm_mip_ui(''Jump'',',...
		'get(get(gcbo,''Parent''),''UserData''),''nrmax'');'],...
	'Interruptible','off','BusyAction','Cancel','Enable',str);
uimenu(h,'Separator','off','Label','goto global maxima',...
	'CallBack',['spm_mip_ui(''Jump'',',...
		'get(get(gcbo,''Parent''),''UserData''),''glmax'');'],...
	'Interruptible','off','BusyAction','Cancel','Enable',str);
uimenu(h,'Separator','on','Label','help',...
	'CallBack','spm_help(''spm_mip_ui'')',...
	'Interruptible','off','BusyAction','Cancel','Enable',str);
set(hMIPim,'UIContextMenu',h)

%-Save handles and data
%-----------------------------------------------------------------------
set(hMIPax,'Tag','hMIPax','UserData',...
	struct(	'hReg',		[],...
		'hMIPxyz',	hMIPxyz,...
		'XYZ',		XYZ,...
		'Z',		Z,...
		'M',		M,...
		'DIM',		DIM,...
		'hXr',		hXr,...
		'hXg',		hXg))

varargout = {hMIPax};



%=======================================================================
case 'getcoords'
%=======================================================================
% xyz = spm_mip_ui('GetCoords',h)
if nargin<2, h=spm_mip_ui('FindMIPax'); else h=varargin{2}; end
varargout = {get(getfield(get(h,'UserData'),'hMIPxyz')	,'UserData')};



%=======================================================================
case 'setcoords'
%=======================================================================
% [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
if nargin<4, hC=0; else, hC=varargin{4}; end
if nargin<3, h=spm_mip_ui('FindMIPax'); else, h=varargin{3}; end
if nargin<2, error('Set co-ords to what!'), else, xyz=varargin{2}; end

MD  = get(h,'UserData');

%-Check validity of coords only when called without a caller handle
%-----------------------------------------------------------------------
if hC<=0
	[xyz,d] = spm_XYZreg('RoundCoords',xyz,MD.M,MD.DIM);
	if d>0 & nargout<2, warning(sprintf(...
	    '%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',...
		mfilename,d)), end
else
	d = [];
end

%-Move marker points, update internal cache in hMIPxyz
%-----------------------------------------------------------------------
spm_mip_ui('PosnMarkerPoints',xyz,h,'r');
set(MD.hMIPxyz,'UserData',reshape(xyz(1:3),3,1))
set(MD.hMIPxyz,'String',{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))})

%-Tell the registry, if we've not been called by the registry...
%-----------------------------------------------------------------------
if (~isempty(MD.hReg) & MD.hReg~=hC), spm_XYZreg('SetCoords',xyz,MD.hReg,h); end

%-Return arguments
%-----------------------------------------------------------------------
varargout = {xyz,d};



%=======================================================================
case 'posnmarkerpoints'
%=======================================================================
% spm_mip_ui('PosnMarkerPoints',xyz,h,r)
if nargin<4, r='r'; else, r=varargin{4}; end
if ~any(strcmp(r,{'r','g'})), error('Invalid pointer colour spec'), end
if nargin<3, h=spm_mip_ui('FindMIPax'); else, h=varargin{3}; end
if nargin<2, xyz = spm_mip_ui('GetCoords',h); else, xyz = varargin{2}; end

%-Get handles of marker points of appropriate colour from UserData of hMIPax
%-----------------------------------------------------------------------
hX = getfield(get(h,'UserData'),['hX',r]);

%-Set marker points
%-----------------------------------------------------------------------
set(hX,'Units','Data')
if length(hX)==1
	set(hX,'Position',[xyz(1), xyz(2), 1])
else
	set(hX(1),'Position',[ Po(1) + xyz(2), Po(2) + xyz(1), 0])
	set(hX(2),'Position',[ Po(1) + xyz(2), Po(3) - xyz(3), 0])
	set(hX(3),'Position',[ Po(4) + xyz(1), Po(3) - xyz(3), 0])
end


%=======================================================================
case 'jump'
%=======================================================================
% [xyz,d] = spm_mip_ui('Jump',h,loc)
if nargin<3, loc='nrvox'; else, loc=varargin{3}; end
if nargin<2, h=spm_mip_ui('FindMIPax'); else, h=varargin{2}; end

%-Get current location & MipData
%-----------------------------------------------------------------------
oxyz = spm_mip_ui('GetCoords',h);
MD   = get(h,'UserData');


%-Compute location to jump to
%-----------------------------------------------------------------------
if isempty(MD.XYZ), loc='dntmv'; end
switch lower(loc), case 'dntmv'
	spm('alert!','No suprathreshold voxels to jump to!',mfilename,0);
	varargout = {oxyz, 0};
	return
case 'nrvox'
	str       = 'nearest suprathreshold voxel';
	[xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,MD.XYZ);
case 'nrmax'
	str       = 'nearest local maxima';
	iM        = inv(MD.M);
	XYZvox    = iM(1:3,:)*[MD.XYZ; ones(1,size(MD.XYZ,2))];
	[null,null,XYZvox,null] = spm_max(MD.Z,XYZvox);
	XYZ       = MD.M(1:3,:)*[XYZvox; ones(1,size(XYZvox,2))];
	[xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,XYZ);
case 'glmax'
	str       = 'global maxima';
	i         = find(MD.Z==max(MD.Z));
	xyz       = MD.XYZ(:,i);
	d         = sqrt(sum((oxyz-xyz).^2));
otherwise
	warning('Unknown jumpmode')
	varargout = {xyz,0};
	return
end

%-Write jump report, jump, and return arguments
%-----------------------------------------------------------------------
fprintf(['\n\t%s:\tJumped %0.2fmm from [%3.0f, %3.0f, %3.0f],\n\t\t\t',...
		'to %s at [%3.0f, %3.0f, %3.0f]\n'],...
	mfilename, d, oxyz, str, xyz)

spm_mip_ui('SetCoords',xyz,h,h);
varargout = {xyz, d};


%=======================================================================
case 'showgreens'
%=======================================================================
% spm_mip_ui('ShowGreens',xyz,h)
if nargin<3, h=get(0,'CurrentFigure'); else, h=varargin{3}; end
h = spm_mip_ui('FindMIPax',h);
if nargin<2, xyz=get(gcbo,'UserData'); else, xyz=varargin{2}; end

spm_mip_ui('PosnMarkerPoints',xyz,h,'g')
set(getfield(get(h,'UserData'),'hXg'),'Visible','on')

%-If called as a callback (xyz from gco), then set WindowButtonUpFcn
%-----------------------------------------------------------------------
if nargin<2
	set(gcbf,'WindowButtonUpFcn',[...
		'set(gcbf,''WindowButtonUpFcn'',''''),',...
		'spm_mip_ui(''HideGreens'',gcbf)'])
end



%=======================================================================
case 'hidegreens'
%=======================================================================
% spm_mip_ui('HideGreens',h)
if nargin<2, h=get(0,'CurrentFigure'); else, h=varargin{2}; end
h = spm_mip_ui('FindMIPax',h);
set(getfield(get(h,'UserData'),'hXg'),'Visible','off')



%=======================================================================
case 'findmipax'
%=======================================================================
% hMIPax = spm_mip_ui('FindMIPax',h)
% Checks / finds hMIPax handles
%-**** h is handle of hMIPax, or figure containing MIP (default gcf)
if nargin<2, h=get(0,'CurrentFigure'); else, h=varargin{2}; end
if isstr(h), h=spm_figure('FindWin',h); end
if ~ishandle(h), error('invalid handle'), end
if ~strcmp(get(h,'Tag'),'hMIPax'), h=findobj(h,'Tag','hMIPax'); end
if isempty(h), error('MIP axes not found'), end
if length(h)>1, error('Multiple MIPs in this figure'), end
varargout = {h};



%=======================================================================
case 'movestart'
%=======================================================================
% spm_mip_ui('MoveStart')
[cO,cF] = gcbo;
hMIPax  = get(cO,'Parent');
MD      = get(hMIPax,'UserData');

%-Store useful quantities in UserData of gcbo, the object to be dragged
%-----------------------------------------------------------------------
set(hMIPax,'Units','Pixels')
set(cO,'UserData',struct(...
	'hReg',		MD.hReg,...
	'xyz',		spm_mip_ui('GetCoords',hMIPax),...
	'MIPaxPos',	get(hMIPax,'Position')*[1,0;0,1;0,0;0,0],...
	'hMIPxyz',	MD.hMIPxyz,...
	'M',		MD.M,...
	'DIM',		MD.DIM,...
	'hX',		MD.hXr))

%-Initiate dragging
%-----------------------------------------------------------------------
if strcmp(get(cF,'SelectionType'),'normal') | isempty(MD.XYZ)
	%-Set Figure callbacks for drop but no drag (DragType 0)
	%---------------------------------------------------------------
	set(MD.hMIPxyz,'Visible','on','String',...
		{'{\bfSPM}{\itmip}','\itPoint & drop...'})
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
		'Interruptible','off')
	set(cF,'Pointer','CrossHair')
elseif strcmp(get(cF,'SelectionType'),'extend')
	%-Set Figure callbacks for drag'n'drop (DragType 1)
	%---------------------------------------------------------------
	set(MD.hMIPxyz,'Visible','on','String',...
		{'{\bfSPM}{\itmip}','\itDynamic drag & drop...'})
	set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',1)',...
		'Interruptible','off')
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'')',...
		'Interruptible','off')
	set(cF,'Pointer','Fleur')
elseif strcmp(get(cF,'SelectionType'),'alt')
	%-Set Figure callbacks for drag'n'drop with co-ord updating (DragType 2)
	%---------------------------------------------------------------
	set(MD.hMIPxyz,'Visible','on','String',...
		{'{\bfSPM}{\itmip}','\itMagnetic drag & drop...'})
	set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',2)',...
		'Interruptible','off')
	set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'')',...
		'Interruptible','off')
	set(cF,'Pointer','Fleur')
end



%=======================================================================
case 'move'
%=======================================================================
% spm_mip_ui('Move',DragType)
if nargin<2, DragType = 2; else, DragType = varargin{2}; end
cF = gcbf;
cO = gco(cF);

%-Get useful data from UserData of gcbo, the object to be dragged
%-----------------------------------------------------------------------
MS  = get(cO,'UserData');
b2d = MS.DIM(3) == 1;

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
if strcmp(sMarker,'hX1r')
	if b2d
		xyz = [d(1); d(2); MS.xyz(3)];
	else
		xyz = [d(2) - Po(2); d(1) - Po(1); MS.xyz(3)];
	end
elseif strcmp(sMarker,'hX2r')
	xyz = [MS.xyz(1); d(1) - Po(1); Po(3) - d(2)];
elseif strcmp(sMarker,'hX3r')
	xyz = [d(1) - Po(4); MS.xyz(2); Po(3) - d(2)];
else
	error('Can''t work out which marker point')
end

%-Round coordinates according to DragType & set in hMIPxyz's UserData
%-----------------------------------------------------------------------
if DragType==0
	xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
elseif DragType==1
	xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
	hMIPax = get(cO,'Parent');
	MD     = get(hMIPax,'UserData');
	i      = spm_XYZreg('FindXYZ',xyz,MD.XYZ);
elseif DragType==2
	hMIPax = get(cO,'Parent');
	MD     = get(hMIPax,'UserData');
	[xyz,i,d] = spm_XYZreg('NearestXYZ',xyz,MD.XYZ);
end
set(MS.hMIPxyz,'UserData',xyz)

%-Move marker points
%-----------------------------------------------------------------------
% spm_mip_ui('PosnMarkerPoints',xyz,hMIPax,'r')
set(MS.hX,'Units','Data')
if b2d
	set(MS.hX,'Position',[xyz(1), xyz(2), 1],'Visible','on')
else
	set(MS.hX(1),'Position',[ Po(1) + xyz(2), Po(2) + xyz(1), 0])
	set(MS.hX(2),'Position',[ Po(1) + xyz(2), Po(3) - xyz(3), 0])
	set(MS.hX(3),'Position',[ Po(4) + xyz(1), Po(3) - xyz(3), 0])
end


%-Update dynamic co-ordinate strings (if appropriate DragType)
%-----------------------------------------------------------------------
if DragType==0
	spm_mip_ui('MoveEnd')
elseif DragType==1
	if isempty(i)
		str = {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz)};
	else
		str = {'{\bfSPM}{\itmip}: ',...
			sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))};
	end
	set(MD.hMIPxyz,'String',str)
elseif DragType==2
	set(MD.hMIPxyz,'String',...
		{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))})
else
	error('Illegal DragType')
end



%=======================================================================
case 'moveend'
%=======================================================================
% spm_mip_ui('MoveEnd')
cF = gcbf;
cO = gco(cF);
hMIPax  = get(cO,'Parent');
MS      = get(cO,'UserData');

%-Reset WindowButton functions, pointer & SPMmip info-string
%-----------------------------------------------------------------------
set(gcbf,'WindowButtonMotionFcn',' ')
set(gcbf,'WindowButtonUpFcn',' ')
set(gcbf,'Pointer','arrow')
set(MS.hMIPxyz,'String',...
	{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',get(MS.hMIPxyz,'UserData'))})

%-Set coordinates after drag'n'drop, tell registry
%-----------------------------------------------------------------------
% don't need to set internal coordinates 'cos 'move' does that
if ~isempty(MS.hReg)
	spm_XYZreg('SetCoords',get(MS.hMIPxyz,'UserData'),MS.hReg,hMIPax);
end



%=======================================================================
otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end
