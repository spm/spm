function varargout = spm_results_ui(varargin)
% User interface for SPM results: Display and analysis of regional effects
% FORMAT [hReg,SPM,VOL,xX,xSDM] = spm_results_ui
%
% hReg   - handle of MIP XYZ registry object
%          (see spm_XYZreg for details)
% SPM    - structure containing SPM, distribution & filtering details
%          (see spm_getSPM.m for contents)
% VOL    - structure containing details of volume analysed
%          (see spm_getSPM.m for contents)
% xX     - Design Matrix structure
%          (see spm_getSPM.m for contents)
% xSDM   - structure containing contents of SPM.mat file
%          ( see spm_spm.m for contents...
%          ( ...except that xX & XYZ are removed
%
% NB: Results section GUI CallBacks use these data structures by name, which
%     therefore *must* be assigned to the correctly named variables.
%_______________________________________________________________________
%
% The SPM results section is for the interactive exploration and
% characterisation of the results of a statistical analysis.
% 
% The user is prompted to select a SPM{Z}, SPM{T} or SPM{F}, that is
% thresholded at user specified levels. The specification of the
% contrasts to use and the height and size thresholds are described in
% spm_getSPM.m. The resulting SPM is then displayed in the graphics
% window as a maximum intensity projection, alongside the design matrix
% and contrasts employed.
%
% The cursors in the MIP can be moved (dragged) to select a particular
% voxel. The three mouse buttons give different drag and drop behaviour:
% Button 1 - point & drop; Button 2 - "dynamic" drag & drop with
% co-ordinate & SPM value updating; Button 3 - "magnetic" drag & drop,
% where the cursor jumps to the nearest suprathreshold voxel in the
% MIP, and shows the value there. (See spm_mip_ui.m, the MIP GUI handling
% function for further details.)
%
% The current voxel specifies the voxel, suprathreshold cluster, or
% orthogonal planes (planes passing through that voxel) for subsequent
% localised utilities.
%
% A control panel in the interactive window enables interactive
% exploration of the results and is divided into volume, cluster, and 
% voxel level utilities:
%
% Volume level:
%	  (i) List - tabulates p values and statistics
%						- see spm_list.m
%	 (ii) 
%	(iii) Write - write out thresholded SPM
%                                    		- see spm_write_filtered.m
%	 (iv) 
%	  (v) Multiplanar - Launches a MultiPlanar lightbox viewing window
%                                               - see spm_????.m
%
% Cluster level:
%         (i) Maxima - tables all local maxima in the suprathreshold region
%             containing the cursor, going onto multiple pages if necessary.
%             The table lists all local maxima, their locations, statistic
%             values and p-values. (Note that the cursor is repositioned to the
%             nearest local maxima by this option.)
%                                  		- see spm_maxima.m
%
%	 (ii) Small Volume Correction (S.V.C) computes a p value corrected
%             for a small specified volume of interest for, and centred on,
%             the current voxel.
%                                 		- see spm_VOI.m
%	(iii) 
%	 (iv) 
%	  (v)
%
% Voxel level:
%         (i) Plot - Graphs of adjusted and fitted activity against various
%             ordinates. (Note that the cursor is repositioned to the nearest
%             voxel with data by this option.) Additionally, writes out 
%             adjusted data to the MatLab command window.
%                                 		- see spm_graph.m
%
%	 (ii) Slices - slices of the thresholded statistic image overlaid
%             on a secondary image chosen by the user. Three transverse
%             slices are shown, being those at the level of the cursor
%             in the z-axis and the two adjacent to it.
%                                 		- see spm_transverse.m
%
%	(iii) Sections - orthogonal sections of the thresholded statistic image
%             overlaid on a secondary image chosen by the user. The sections
%             are through the cursor position.
%                                    		- see spm_sections.m
%	 (iv) 
%	  (v) 
%
% Graphics appear in the bottom half of the graphics window, additional
% controls and questions appearing in the interactive window.
%
%                           ----------------
%
% The MIP uses a template outline in Talairach space. Consequently for
% the results section to display properly the input images to the
% statistics section should either be in Talairach space (with the
% ORIGIN correctly specified), or the ORIGIN header fields should be
% set to the voxel coordinates of the anterior commissure in the input
% images. See spm_format.man ("Data Format" in the help facility) for
% further details of the Analyze image format used by SPM.
%
% Similarly, secondary images should be aligned with the input images
% used for the statistical analysis. In particular the ORIGIN must
% correspond to (0,0,0) in XYZ, the vector of locations.
%
%                           ----------------
%
% In addition to setting up the results section, spm_results_ui.m sets
% up the results section GUI and services the CallBacks. FORMAT
% specifications for embedded CallBack functions are given in the main
% body of the code.
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%
SCCSid = '%I%';

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `spm_figure Create` is           )
%( equivalent to `spm('Create')`.                                      )
%
% spm_results_ui sets up and handles the SPM results graphical user
% unterface, initialising an XYZ registry (see spm_XYZreg.m) to co-ordinate
% locations between various location controls.
%
%_______________________________________________________________________
%
% FORMAT hReg = spm_results_ui('SetupGUI',M,DIM,Finter)
% Setup results GUI in Interactive window
% M       - 4x4 transformation matrix relating voxel to "real" co-ordinates
% DIM     - 3 vector of image X, Y & Z dimensions
% Finter  - handle (or 'Tag') of Interactive window (default 'Interactive')
% hReg    - handle of XYZ registry object
%
% FORMAT [hFvolB,hFclusB,hFvoxB] = ...
%	spm_results_ui('DrawButts',hReg,DIM,Finter,WS,FS)
% Draw GUI buttons
% hReg    - handle of XYZ registry object
% DIM     - 3 vector of image X, Y & Z dimensions
% Finter  - handle of Interactive window
% WS      - WinScale  [Default spm('WinScale') ]
% FS      - FontSizes [Default spm('FontSizes')]
% hFvolB  - handle to frame enclosing volume level buttons
% hFclusB - handle to frame enclosing cluster level buttons
% hFvoxB  - handle to frame enclosing voxel level buttons
%
% FORMAT hFxyz = spm_results_ui('DrawXYZgui',M,DIM,xyz,Finter)
% Setup editable XYZ control widgets at foot of Interactive window
% M      - 4x4 transformation matrix relating voxel to "real" co-ordinates
% DIM    - 3 vector of image X, Y & Z dimensions
% xyz    - Initial xyz location {mm}
% Finter - handle of Interactive window
% hFxyz  - handle of XYZ control - the frame containing the edit widgets
%
% FORMAT spm_results_ui('EdWidCB')
% Callback for editable XYZ control widgets
%
% FORMAT xyz = spm_results_ui('GetCoords',hFxyz)
% Get current co-ordinates from editable XYZ control
% hFxyz - handle of frame enclosing widgets - the Tag object for this control
% xyz   - current co-ordinates {mm}
% NB: When using the results section, should use XYZregistry to get/set location
%
% FORMAT [xyz,d] = spm_results_ui('SetCoords',xyz,hFxyz,hC)
% Set co-ordinates to XYZ widget
% xyz   - (Input) desired co-ordinates {mm}
% hFxyz - handle of XYZ control - the frame containing the edit widgets
% hC    - handle of calling object, if used as a callback. [Default 0]
% xyz   - (Output) Desired co-ordinates are rounded to nearest voxel if hC
%         is not specified, or is zero. Otherwise, caller is assummed to
%         have checked verity of desired xyz co-ordinates. Output xyz returns
%         co-ordinates actually set {mm}.
% d     - Euclidean distance between desired and set co-ordinates.
% NB: When using the results section, should use XYZregistry to get/set location
%
% FORMAT hFxyz = spm_results_ui('FindXYZframe',h)
% Find/check XYZ edit widgets frame handle, 'Tag'ged 'hFxyz'
% h     - handle of frame enclosing widgets, or containing figure [default gcf]
%         If isstr(h), then uses spm_figure('FindWin',h) to locate named figures
% hFxyz - handle of confirmed XYZ editable widgets control
%         Errors if hFxyz is not an XYZ widget control, or a figure containing
%         a unique such control
%
% FORMAT spm_results_ui('PlotUi',hAx)
% GUI for adjusting plot attributes - Sets up controls just above results GUI
% hAx - handle of axes to work with
%
% FORMAT spm_results_ui('PlotUiCB')
% CallBack handler for Plot attribute GUI
%
% FORMAT Fgraph = spm_results_ui('Clear',F,RNP)
% Clears results subpane of Graphics window, deletes PageControls, sets new axes
% F      - handle of Graphics window [Default spm_figure('FindWin','Graphics')]
% RNP    - If specified, respects 'NextPlot' axes attribute, and won't delete
%          axes with NextPlot set to 'add' (which is set by `hold on`)
% Fgraph - handle of Graphics window
%
% FORMAT hMP = spm_results_ui('LaunchMP',M,DIM,hReg,hBmp)
% Prototype callback handler for integrating MultiPlanar toolbox
%
% FORMAT spm_results_ui('Delete',h)
% deletes HandleGraphics objects, but only if they're valid, thus avoiding
% warning statements from MatLab!
%_______________________________________________________________________


%-Condition arguments
%-----------------------------------------------------------------------
if nargin==0, Action='SetUp'; else, Action=varargin{1}; end


%=======================================================================
switch lower(Action), case 'setup'                      %-Set up results
%=======================================================================
%-Initialise 
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
FS    = spm('FontSizes');

%-Get thresholded SPM data and parameters of design
%=======================================================================
[SPM,VOL,xX,xCon,xSDM] = spm_getSPM;


%-Setup Results User Interface; Display MIP, design matrix & parameters
%=======================================================================
spm('FigName',['SPM{',SPM.STAT,'}: Results'],Finter,CmdLine);


%-Setup results GUI
%-----------------------------------------------------------------------
hReg   = spm_results_ui('SetupGUI',VOL.M,VOL.DIM,Finter);


%-Setup Maximium intensity projection (MIP) & register
%-----------------------------------------------------------------------
hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
hMIPax = spm_mip_ui(SPM.Z,SPM.XYZmm,VOL.M,VOL.DIM,hMIPax);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
text(260,260,['SPM\{',SPM.STATstr,'\}'],...
	'Interpreter','TeX',...
	'FontSize',FS(14),'Fontweight','Bold',...
	'Parent',hMIPax)


%-Print comparison title
%-----------------------------------------------------------------------
hTitAx = axes('Parent',Fgraph,...
		'Position',[0.02 0.95 0.96 0.02],...
		'Visible','off');

text(0.5,0,SPM.title,'Parent',hTitAx,...
	'HorizontalAlignment','center',...
	'VerticalAlignment','baseline',...
	'FontWeight','Bold','FontSize',FS(14))


%-Print SPMresults: Results directory & thresholding info
%-----------------------------------------------------------------------
hResAx = axes('Parent',Fgraph,...
		'Position',[0.05 0.55 0.45 0.05],...
		'DefaultTextVerticalAlignment','baseline',...
		'DefaultTextFontSize',FS(9),...
		'DefaultTextColor',[1,1,1]*.7,...
		'Units','points',...
		'Visible','off');
AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
h = text(0,24,'SPMresults:','Parent',hResAx,...
	'FontWeight','Bold','FontSize',FS(14));
text(get(h,'Extent')*[0;0;1;0],24,spm_str_manip(SPM.swd,'a30'),'Parent',hResAx)
text(0,12,sprintf('Height threshold %c = %0.2f',SPM.STAT,SPM.u),'Parent',hResAx)
text(0,0,sprintf('Extent threshold k = %0.0f resels',SPM.k),'Parent',hResAx)


%-Plot design matrix
%-----------------------------------------------------------------------
hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
hDesMtxIm = image((xX.nKX+1)*32);
set(hDesMtx,'YTick',spm_DesRep('ScanTick',size(xX.nKX,1),24),'TickDir','out')
xlabel('Design matrix')
set(hDesMtxIm,'UserData',...
	struct('X',xX.xKXs.X,'fnames',{{xSDM.VY.fname}'},'Xnames',{xX.Xnames}))
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')


%-Plot contrasts
%-----------------------------------------------------------------------
nPar = size(xX.X,2);
xx   = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
nCon = length(SPM.Ic);

dy    = 0.15/max(nCon,2);
for ii = nCon:-1:1
    axes('Position',[0.65 (0.80 + dy*(nCon-ii+.1)) 0.25 dy*.9])
    if xCon(SPM.Ic(ii)).STAT == 'T' & size(xCon(SPM.Ic(ii)).c,2)==1
	%-Single vector contrast for SPM{t} - bar
	yy = [zeros(1,nPar);repmat(xCon(SPM.Ic(ii)).c',2,1);zeros(1,nPar)];
	patch(xx,yy,[1,1,1]*.5)
	set(gca,'Box','off','XTick',[],'YTick',[])
	set(gca,'Tag','ConGrphAx')
    else
	%-F-contrast - image
	sca = 1/max(abs(xCon(SPM.Ic(ii)).c(:)));
	image((xCon(SPM.Ic(ii)).c'*sca+1)*32)
	set(gca,'Box','on','XTick',[],'YTick',[])
	set(gca,'Tag','ConGrphAx')
    end
    ylabel(num2str(SPM.Ic(ii)))
    axis tight
end
title('Contrast(s)')


%-Store handles of results section Graphics window objects
%-----------------------------------------------------------------------
H  = get(Fgraph,'Children');
H  = findobj(H,'flat','HandleVisibility','on');
H  = findobj(H);
Hv = get(H,'Visible');
set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))


%-Finished results setup
%-----------------------------------------------------------------------
varargout = {hReg,SPM,VOL,xX,xCon,xSDM};
spm('Pointer','Arrow')



%=======================================================================
case 'setupgui'                             %-Set up results section GUI
%=======================================================================
% hReg = spm_results_ui('SetupGUI',M,DIM,Finter)
if nargin<4, Finter='Interactive'; else, Finter=varargin{4}; end
if nargin<3, error('Insufficient arguments'), end
M      = varargin{2};
DIM    = varargin{3};

Finter = spm_figure('GetWin',Finter);
WS     = spm('WinScale');
FS     = spm('FontSizes');

%-Create frame for Results GUI objects; "Help"; "Exit", & "Clear" buttons...
%-----------------------------------------------------------------------
hReg    = uicontrol(Finter,'Style','Frame','BackgroundColor',spm('Colour'),...
		'Position',[001 001 400 175].*WS);
hFResUi = uicontrol(Finter,'Style','Frame','Position',[008 007 387 163].*WS);

hClear = uicontrol(Finter,'Style','PushButton','String','clear',...
	'ToolTipString','clears lower results subpane of Graphics window',...
	'FontSize',FS(10),'ForegroundColor','b',...
	'Callback',[...
		'spm_results_ui(''Clear''); ',...
		'spm_input(''!DeleteInputObj'')'],...
	'Interruptible','on','Enable','on',...
	'DeleteFcn','clc,spm_clf(''Graphics'')',...
	'Position',[325 028 045 018].*WS);
hExit  = uicontrol(Finter,'Style','PushButton','String','exit',...
	'ToolTipString','exit the results section',...
	'FontSize',FS(10),'ForegroundColor','r',...
	'Callback','spm_clf(''Interactive''), spm_clf(''Graphics''), clear',...
	'Interruptible','on','Enable','on',...
	'Position',[325 010 045 018].*WS);
hHelp  = uicontrol(Finter,'Style','PushButton','String','?',...
	'ToolTipString','results section help',...
	'FontSize',FS(10),'ForegroundColor','g',...
	'Callback','spm_help(''spm_results_ui'')',...
	'Interruptible','on','Enable','on',...
	'Position',[370 010 020 036].*WS);


%-Initialise registry in hReg frame object
%-----------------------------------------------------------------------
[hReg,xyz] = spm_XYZreg('InitReg',hReg,M,DIM,[0;0;0]);

%-Setup editable XYZ widgets & cross register with registry
%-----------------------------------------------------------------------
hFxyz      = spm_results_ui('DrawXYZgui',M,DIM,xyz,Finter);
spm_XYZreg('XReg',hReg,hFxyz,'spm_results_ui');

%-Set up buttons for Volume, Cluster & Voxel level utility functions
%-----------------------------------------------------------------------
[hFvolB,hFclusB,hFvoxB] = spm_results_ui('DrawButts',hReg,DIM,Finter,WS,FS);

varargout  = {hReg};



%=======================================================================
case 'drawbutts'    %-Draw results section buttons in Interactive window
%=======================================================================
% [hFvolB,hFclusB,hFvoxB] = spm_results_ui('DrawButts',hReg,DIM,Finter,WS,FS)
if nargin<3, error('Insufficient arguments'), end
hReg = varargin{2};
DIM  = varargin{3};
if nargin<4,  Finter=spm_figure('FindWin','Interactive');
	else, Finter=varargin{4}; end
if nargin < 5, WS = spm('WinScale');	else,	WS = varargin{5}; end
if nargin < 6, FS = spm('FontSizes');	else,	FS = varargin{6}; end

if DIM(3) > 1, s3D = 'on'; else, s3D = 'off'; end


%-Volume level functions
%-----------------------------------------------------------------------
hFvolB = uicontrol(Finter,'Style','Frame','Position',[010 130 380 030].*WS);
uicontrol(Finter,'Style','Text','String','Volume level functions',...
	'Position',[020 152 125 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','p values','FontSize',FS(10),...
	'ToolTipString','tabulate summary of p-values & statistics',...
	'Callback','spm_list(SPM,VOL,[],[],''volume summary...'')',...
	'Interruptible','on','Enable','on',...
	'Position',[015 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[090 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','write','FontSize',FS(10),...
	'ToolTipString','Write current current filtered SPM',...
	'Callback',['spm_write_filtered(SPM.Z,SPM.XYZ,VOL.DIM,VOL.M,',...
		'sprintf(''SPM{%c}-filtered: u = %5.3f, k = %d'',',...
			'SPM.STAT,SPM.u,SPM.k));'],...
	'Interruptible','on','Enable','on',...
	'Position',[165 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','multiplane','FontSize',FS(10),...
	'ToolTipString','Launch MultiPlanar viewing toolbox',...
	'Tag','hB_MultiPlanar',...
	'Callback',...
		'hMP = spm_results_ui(''LaunchMP'',VOL.M,VOL.DIM,hReg,gcbo);',...
	'Interruptible','on','Enable','off',...
	'Position',[315 135 070 020].*WS)

%-Cluster level functions
%-----------------------------------------------------------------------
hFclusB = uicontrol(Finter,'Style','Frame','Position',[010 090 380 030].*WS);
uicontrol(Finter,'Style','Text','String','Cluster level functions',...
	'Position',[020 112 125 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','maxima','FontSize',FS(10),...
	'ToolTipString',['tabulate p-values & statistics for local maxima ',...
		'of current cluster'],...
	'Callback','spm_maxima(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[015 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','S.V.C.','FontSize',FS(10),...
	'ToolTipString',['Small Volume Correction - corrected p-values ',...
		'for a small search region centered on this voxel'],...
	'Callback','spm_VOI(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[090 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[165 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[315 095 070 020].*WS)

%-Voxel level functions
%-----------------------------------------------------------------------
hFvoxB = uicontrol(Finter,'Style','Frame','Position',[010 050 380 030].*WS);
uicontrol(Finter,'Style','Text','String','Voxel level functions',...
	'Position',[020 072 115 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','plot','FontSize',FS(10),...
	'ToolTipString','plotting of parameters & data at this voxel',...
	'Callback','[Y,y,beta,SE] = spm_graph(SPM,VOL,xX,xCon,xSDM,hReg);',...
	'Interruptible','on','Enable','on',...
	'Position',[015 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','slices','FontSize',FS(10),...
	'ToolTipString','overlay SPM on three slices of another image',...
	'Callback','spm_transverse(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[090 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','sections','FontSize',FS(10),...
	'ToolTipString',['overlay filtered SPM on orthogonal sections of',...
		'another image'],...
	'Callback','spm_sections(SPM,VOL,hReg)',...
	'Interruptible','on','Enable',s3D,...
	'Position',[165 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',FS(10),...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[315 055 070 020].*WS)

varargout = {hFvolB,hFclusB,hFvoxB};

%=======================================================================
case 'drawxyzgui'                                    %-Draw XYZ GUI area
%=======================================================================
% hFxyz = spm_results_ui('DrawXYZgui',M,DIM,xyz,Finter)
if nargin<5,  Finter=spm_figure('FindWin','Interactive');
	else, Finter=varargin{5}; end
if nargin<4, xyz=[0;0;0]; else, xyz=varargin{4}; end
if nargin<3, error('Insufficient arguments'), end
DIM     = varargin{3};
M       = varargin{2};

xyz     = spm_XYZreg('RoundCoords',xyz,M,DIM);

%-Locate windows etc...
%-----------------------------------------------------------------------
WS      = spm('WinScale');
FS      = spm('FontSizes');

%-Create & link control objects
%-----------------------------------------------------------------------
hFxyz = uicontrol(Finter,'Style','Frame','Position',[010 010 310 030].*WS);
uicontrol(Finter,'Style','Text','String','Co-ordinates',...
	'Position',[020 032 080 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')


uicontrol(Finter,'Style','Text','String','x =',...
	'Position',[020 015 020 018].*WS,...
	'FontName','Times','FontSize',FS(10),...
	'HorizontalAlignment','Center');
hX   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(1)),...
	'ToolTipString','enter x-coordinate',...
	'Position',[040 015 070 020].*WS,...
	'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hX',...
	'Callback','spm_results_ui(''EdWidCB'')');

uicontrol(Finter,'Style','Text','String','y =',...
	'Position',[120 015 020 018].*WS,...
	'FontName','Times','FontSize',FS(10),...
	'HorizontalAlignment','Center')
hY   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(2)),...
	'ToolTipString','enter y-coordinate',...
	'Position',[140 015 070 020].*WS,...
	'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hY',...
	'Callback','spm_results_ui(''EdWidCB'')');

uicontrol(Finter,'Style','Text','String','z =',...
	'Position',[220 015 020 018].*WS,...
	'FontName','Times','FontSize',FS(10),...
	'HorizontalAlignment','Center')
hZ   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(3)),...
	'ToolTipString','enter z-coordinate',...
	'Position',[240 015 070 020].*WS,...
	'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hZ',...
	'Callback','spm_results_ui(''EdWidCB'')');

set(hFxyz,'Tag','hFxyz','UserData',struct(...
	'hReg',	[],...
	'M',	M,...
	'DIM',	DIM,...
	'hX',	hX,...
	'hY',	hY,...
	'hZ',	hZ,...
	'xyz',	xyz	));

set([hX,hY,hZ],'UserData',hFxyz)

%uicontrol(Finter,'Style','Text','String',spm('Ver'),...
%	'Position',[310 005 080 020].*WS,...
%	'FontName','Times','FontWeight','Bold',...
%	'HorizontalAlignment','Center',...
%	'ForegroundColor',[1 1 1]*.6)

varargout = {hFxyz};



%=======================================================================
case 'edwidcb'                           %-Callback for editable widgets
%=======================================================================
% spm_results_ui('EdWidCB')

hC   = gcbo;
d    = find(strcmp(get(hC,'Tag'),{'hX','hY','hZ'}));
hF   = get(hC,'UserData');
UD   = get(hF,'UserData');
xyz  = UD.xyz;
nxyz = xyz;

o = evalin('base',['[',get(hC,'String'),']'],'sprintf(''error'')');
if ischar(o) | length(o)>1
	warning(sprintf('%s: Error evaluating ordinate:\n\t%s',...
		mfilename,lasterr))
else
	nxyz(d) = o;
	nxyz = spm_XYZreg('RoundCoords',nxyz,UD.M,UD.DIM);
end

if abs(xyz(d)-nxyz(d))>0
	UD.xyz = nxyz; set(hF,'UserData',UD)
	if ~isempty(UD.hReg), spm_XYZreg('SetCoords',nxyz,UD.hReg,hF); end
end

set(hC,'String',sprintf('%.3f',nxyz(d)))



%=======================================================================
case 'getcoords'              % Get current co-ordinates from XYZ widget
%=======================================================================
% xyz = spm_results_ui('GetCoords',hFxyz)
if nargin<2, hFxyz='Interactive'; else, hFxyz=varargin{2}; end
hFxyz     = spm_results_ui('FindXYZframe',hFxyz);
varargout = {getfield(get(hFxyz,'UserData'),'xyz')};



%=======================================================================
case 'setcoords'                        % Set co-ordinates to XYZ widget
%=======================================================================
% [xyz,d] = spm_results_ui('SetCoords',xyz,hFxyz,hC)
if nargin<4, hC=0; else, hC=varargin{4}; end
if nargin<3, hFxyz=spm_results_ui('FindXYZframe'); else, hFxyz=varargin{3}; end
if nargin<2, error('Set co-ords to what!'), else, xyz=varargin{2}; end

%-If this is an internal call, then don't do anything
if hFxyz==hC, return, end

UD = get(hFxyz,'UserData');

%-Check validity of coords only when called without a caller handle
%-----------------------------------------------------------------------
if hC<=0
	[xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.DIM);
	if d>0 & nargout<2, warning(sprintf(...
	    '%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',...
		mfilename,d)), end
else
	d = [];
end

%-Update xyz information & widget strings
%-----------------------------------------------------------------------
UD.xyz = xyz; set(hFxyz,'UserData',UD)
set(UD.hX,'String',sprintf('%.3f',xyz(1)))
set(UD.hY,'String',sprintf('%.3f',xyz(2)))
set(UD.hZ,'String',sprintf('%.3f',xyz(3)))

%-Tell the registry, if we've not been called by the registry...
%-----------------------------------------------------------------------
if (~isempty(UD.hReg) & UD.hReg~=hC)
	spm_XYZreg('SetCoords',xyz,UD.hReg,hFxyz);
end

%-Return arguments
%-----------------------------------------------------------------------
varargout = {xyz,d};



%=======================================================================
case 'findxyzframe'                                   % Find hFxyz frame
%=======================================================================
% hFxyz = spm_results_ui('FindXYZframe',h)
% Sorts out hFxyz handles
if nargin<2, h='Interactive'; else, h=varargin{2}; end
if isstr(h), h=spm_figure('FindWin',h); end
if ~ishandle(h), error('invalid handle'), end
if ~strcmp(get(h,'Tag'),'hFxyz'), h=findobj(h,'Tag','hFxyz'); end
if isempty(h), error('XYZ frame not found'), end
if length(h)>1, error('Multiple XYZ frames found'), end
varargout = {h};



%=======================================================================
case 'plotui'                                %-GUI for plot manipulation
%=======================================================================
% spm_results_ui('PlotUi',hAx)
if nargin<2, hAx=gca; else, hAx=varargin{2}; end

WS = spm('WinScale');
FS = spm('FontSizes');
Finter=spm_figure('FindWin','Interactive');
figure(Finter)

%-Check there aren't already controls!
%-----------------------------------------------------------------------
hGraphUI = findobj(Finter,'Tag','hGraphUI');
if ~isempty(hGraphUI)			%-Controls exist
	hBs = get(hGraphUI,'UserData');
	if hAx==get(hBs(1),'UserData')	%-Controls linked to these axes
		return
	else				%-Old controls remain
		delete(findobj(Finter,'Tag','hGraphUIbg'))
	end
end

%-Frames & text
%-----------------------------------------------------------------------
hGraphUIbg = uicontrol(Finter,'Style','Frame','Tag','hGraphUIbg',...
		'BackgroundColor',spm('Colour'),...
		'Position',[001 181 400 055].*WS);
hGraphUI = uicontrol(Finter,'Style','Frame','Tag','hGraphUI',...
		'Position',[008 187 387 043].*WS);
hGraphUIButtsF = uicontrol(Finter,'Style','Frame',...
		'Position',[010 190 380 030].*WS);
hText = uicontrol(Finter,'Style','Text','String','Plot controls',...
	'Position',[020 212 070 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w');

%-Controls
%-----------------------------------------------------------------------
h1 = uicontrol(Finter,'Style','CheckBox','String','Hold',...
	'ToolTipString','toggle hold to overlay plots',...
	'FontSize',FS(10),...
	'Value',strcmp(get(hAx,'NextPlot'),'add'),...
	'Callback',[...
		'if get(gcbo,''Value''), ',...
		    'set(get(gcbo,''UserData''),''NextPlot'',''add''), ',...
		'else, ',...
		    'set(get(gcbo,''UserData''),''NextPlot'',''replace''), ',...
		'end'],...
	'Interruptible','on','Enable','on',...
	'Position',[015 195 070 020].*WS);
h2 = uicontrol(Finter,'Style','CheckBox','String','Grid',...
	'ToolTipString','toggle axes grid',...
	'FontSize',FS(10),...
	'Value',strcmp(get(hAx,'XGrid'),'on'),...
	'Callback',[...
		'if get(gcbo,''Value''), ',...
			'set(get(gcbo,''UserData''),''XGrid'',''on'','...
		    		'''YGrid'',''on'',''ZGrid'',''on''), ',...
		'else, ',...
			'set(get(gcbo,''UserData''),''XGrid'',''off'','...
		    		'''YGrid'',''off'',''ZGrid'',''off''), ',...
		'end'],...
	'Interruptible','on','Enable','on',...
	'Position',[090 195 070 020].*WS);
h3 = uicontrol(Finter,'Style','CheckBox','String','Box',...
	'ToolTipString','toggle axes box',...
	'FontSize',FS(10),...
	'Value',strcmp(get(hAx,'Box'),'on'),...
	'Callback',[...
		'if get(gcbo,''Value''), ',...
		    'set(get(gcbo,''UserData''),''Box'',''on''), ',...
		'else, ',...
		    'set(get(gcbo,''UserData''),''Box'',''off''), ',...
		'end'],...
	'Interruptible','on','Enable','on',...
	'Position',[165 195 070 020].*WS);
h4 = uicontrol(Finter,'Style','PopUp',...
	'ToolTipString','edit axis text annotations',...
	'FontSize',FS(10),...
	'String','text|Title|Xlabel|Ylabel',...
	'Callback','spm_results_ui(''PlotUiCB'')',...
	'Interruptible','on','Enable','on',...
	'Position',[240 195 070 020].*WS);
h5 = uicontrol(Finter,'Style','PopUp',...
	'ToolTipString','change various axes attributes',...
	'FontSize',FS(10),...
	'String','attrib|LineWidth|XLim|YLim|handle',...
	'Callback','spm_results_ui(''PlotUiCB'')',...
	'Interruptible','off','Enable','on',...
	'Position',[315 195 070 020].*WS);

%-Handle storage for linking, and DeleteFcns for linked deletion
%-----------------------------------------------------------------------
set(hGraphUI,'UserData',[h1,h2,h3,h4,h5])
set([h1,h2,h3,h4,h5],'UserData',hAx)

set(hGraphUIbg,'UserData',...
	[hGraphUI,hGraphUIButtsF,hText,h1,h2,h3,h4,h5],...
	'DeleteFcn','spm_results_ui(''Delete'',get(gcbo,''UserData''))')
set(hAx,'UserData',hGraphUIbg,...
	'DeleteFcn','spm_results_ui(''Delete'',get(gcbo,''UserData''))')




%=======================================================================
case 'plotuicb'
%=======================================================================
% spm_results_ui('PlotUiCB')
hPM = gcbo;
v   = get(hPM,'Value');
if v==1, return, end
str = cellstr(get(hPM,'String'));
str = str{v};

hAx = get(hPM,'UserData');
switch str
case 'Title'
	h = get(hAx,'Title');
	set(h,'String',spm_input('Enter title:',-1,'s',get(h,'String')))
case 'Xlabel'
	h = get(hAx,'Xlabel');
	set(h,'String',spm_input('Enter X axis label:',-1,'s',get(h,'String')))
case 'Ylabel'
	h = get(hAx,'Ylabel');
	set(h,'String',spm_input('Enter Y axis label:',-1,'s',get(h,'String')))
case 'LineWidth'
	lw = spm_input('Enter LineWidth',-1,'e',get(hAx,'LineWidth'));
	set(hAx,'LineWidth',lw)
case 'XLim'
	XLim = spm_input('Enter XLim',-1,'e',get(hAx,'XLim'));
	set(hAx,'XLim',XLim)
case 'YLim'
	YLim = spm_input('Enter YLim',-1,'e',get(hAx,'YLim'));
	set(hAx,'YLim',YLim)
case 'handle'
	varargout={hAx};
otherwise
	warning(['Unknown action: ',str])
end

set(hPM,'Value',1)


%=======================================================================
case {'clear','clearpane'}                       %-Clear results subpane
%=======================================================================
% Fgraph = spm_results_ui('Clear',F,mode)
% mode 1 [default] usual, mode 0 - clear & hide Res stuff, 2 - RNP
if strcmp(lower(Action),'clearpane')
	warning('''ClearPane'' action is grandfathered, use ''Clear'' instead')
end

if nargin<3, mode=1; else, mode=varargin{3}; end
if nargin<2, F='Graphics'; else, F=varargin{2}; end
F = spm_figure('FindWin',F);

%-Clear input objects from 'Interactive' window
%-----------------------------------------------------------------------
%spm_input('!DeleteInputObj')


%-Get handles of objects in Graphics window & note permanent results objects
%-----------------------------------------------------------------------
H = get(F,'Children');				%-Get contents of window
H = findobj(H,'flat','HandleVisibility','on');	%-Lose GUI components
h = findobj(H,'flat','Tag','PermRes');		%-Look for 'Res' obj
if ~isempty(h)					% (which has perm. res. obj. h's)
	tmp  = get(h,'UserData');
	HR   = tmp.H;
	HRv  = tmp.Hv;
elseif mode==0
	HR   = [];
	HRv  = {};
else
	%-Can't find stored handles of permanent results objects
	% (Work out as everything in top half of figure)
	%-This won't work properly with "hiding" of permanent results objects
	warning('Can''t find cache of permanent results objects - deriving')
	un   = get(H,'Units');
	set(H,'Units','normalized')
	tmp  = get(H,'Position');
	if length(H)==1
		un = {un};
		if tmp(2)>0.5
			HR=H;  HRv = {get(HR,'Visible')};
		else
			HR=[]; HRv = {};
		end
	elseif ~isempty(tmp)
		tmp  = cat(1,tmp{:});
		HR   = H(tmp(:,2)>0.5);
		HRv  = get(HR,'Visible');
	else
		HR   = [];
		HRv  = {};
	end
	set(H,{'Units'},un)
end
H = setdiff(H,HR);				%-Drop permanent results obj



%-Delete stuff as appropriate
%-----------------------------------------------------------------------
if mode==2	%-Don't delete axes with NextPlot 'add'
	H = setdiff(H,findobj(H,'flat','Type','axes','NextPlot','add'));
end

delete(H)

if mode==0	%-Hide the permanent results section stuff
	set(HR,'Visible','off')
else
	set(HR,{'Visible'},HRv)
end


%=======================================================================
case 'launchmp'                             %-Launch multiplanar toolbox
%=======================================================================
% hMP = spm_results_ui('LaunchMP',M,DIM,hReg,hBmp)
if nargin<5, hBmp = gcbo; else, hBmp = varargin{5}; end
hReg = varargin{4};
DIM  = varargin{3};
M    = varargin{2};

%-Check for existing MultiPlanar toolbox
hMP = get(hBmp,'UserData');
if ishandle(hMP)
	figure(spm_figure('ParentFig',hMP))
	varargout = {hMP};
	return
end

%-Initialise and cross-register MultiPlanar toolbox
hMP = spm_XYZreg_Ex2('Create',M,DIM);
spm_XYZreg('Xreg',hReg,hMP,'spm_XYZreg_Ex2');

%-Setup automatic deletion of MultiPlanar on deletion of results controls
set(hBmp,'Enable','on','UserData',hMP)
set(hBmp,'DeleteFcn','spm_results_ui(''delete'',get(gcbo,''UserData''))')

varargout = {hMP};



%=======================================================================
case 'delete'                            %-Delete HandleGraphics objects
%=======================================================================
% spm_results_ui('Delete',h)
h = varargin{2};
delete(h(ishandle(h)));


%=======================================================================
otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end

