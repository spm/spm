function varargout = spm_results_ui(varargin)
% User interface for SPM results section
% FORMAT varargout = spm_results_ui(varargin)
%_______________________________________________________________________
%
% spm_results_ui sets up and handles the SPM results graphical user
% unterface, initialising an XYZ registry (see spm_XYZreg.m) to co-ordinate
% locations between various location controls
%
% See spm_results.m for help on using the results section GUI...
% Brief programmers help & specifications are given below...
%_______________________________________________________________________
%
% FORMAT hReg = spm_results_ui('SetupGUI',M,D)
% Setup results GUI in Interactive window
% M    - 4x4 transformation matrix relating voxel to "real" co-ordinates
% D    - 3 vector of image X, Y & Z dimensions (DIM)
% hReg - Handle of XYZ registry object
%
% FORMAT [hFvolB,hFclusB,hFvoxB] = ...
%	spm_results_ui('DrawButts',hReg,D,Finter,WS,Fs)
% Draw GUI buttons
% hReg    - Handle of XYZ registry object
% D       - 3 vector of image X, Y & Z dimensions (DIM)
% Finter  - handle of Interactive window
% WS      - WinScale [Default spm('GetWinScale')]
% Fs      - FontSize [Default 2*round(10*min(WS)/2)]
% hFvolB  - handle to frame enclosing volume level buttons
% hFclusB - handle to frame enclosing cluster level buttons
% hFvoxB  - handle to frame enclosing voxel level buttons
%
% FORMAT hFxyz = spm_results_ui('DrawXYZgui',M,D,xyz,Finter)
% Setup editable XYZ control widgets at foot of Interactive window
% M      - 4x4 transformation matrix relating voxel to "real" co-ordinates
% D      - 3 vector of image X, Y & Z dimensions (DIM)
% xyz    - Initial xyz location
% Finter - handle of Interactive window
% hFxyz  - handle of XYZ control - the frame containing the edit widgets
%
% FORMAT spm_results_ui('EdWidCB')
% Callback for editable XYZ control widgets
%
% FORMAT xyz = spm_results_ui('GetCoords',hFxyz)
% Get current co-ordinates from editable XYZ control
% hFxyz - handle of frame enclosing widgets - the Tag object for this control
% xyz   - current co-ordinates
% NB: When using the results section, should use XYZregistry to get/set location
%
% FORMAT [xyz,d] = spm_results_ui('SetCoords',xyz,hFxyz,hC)
% Set co-ordinates to XYZ widget
% xyz   - (Input) desired co-ordinates
% hFxyz - handle of XYZ control - the frame containing the edit widgets
% hC    - Handle of calling object, if used as a callback. [Default 0]
% xyz   - (Output) Desired co-ordinates are rounded to nearest voxel if hC
%         is not specified, or is zero. Otherwise, caller is assummed to
%         have checked verity of desired xyz co-ordinates. Output xyz returns
%         co-ordinates actually set.
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
% FORMAT Fgraph = spm_results_ui('ClearPane',F,RNP)
% Clears results subpane of Graphics window, deletes PageControls, sets new axes
% F      - handle of Graphics window [Default spm_figure('FindWin','Graphics')]
% RNP    - If specified, respects 'NextPlot' axes attribute, and won't delete
%          axes with NextPlot set to 'add' (which is set by `hold on`)
% Fgraph - Handle of Graphics window
%
% FORMAT hMP = spm_results_ui('LaunchMP',M,DIM,hReg,hBmp)
% Prototype callback handler for integrating MultiPlanar toolbox
%
% FORMAT spm_results_ui('Delete',h)
% deletes HandleGraphics objects, but only if they're valid, thus avoiding
% warning statements from MatLab!
%_______________________________________________________________________
% %W% Andrew Holmes %E%




%=======================================================================
switch lower(varargin{1}), case 'setupgui'	%-Set up results GUI
%=======================================================================
% hReg = spm_results_ui('SetupGUI',M,D)
if nargin<3, error('Insufficient arguments'), end
M      = varargin{2};
D      = varargin{3};

Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), error('No Interactive window'), end
set(Finter,'Name','SPM results')
WS     = spm('GetWinScale');
Fs     = 2*round(10*min(WS)/2);

%-Create frame for Results GUI objects; "Help"; "Exit", & "Clear" buttons...
%-----------------------------------------------------------------------
hReg    = uicontrol(Finter,'Style','Frame','BackgroundColor',spm('Colour'),...
		'Position',[001 001 400 175].*WS);
hFResUi = uicontrol(Finter,'Style','Frame','Position',[008 007 387 163].*WS);

hClear = uicontrol(Finter,'Style','PushButton','String','clear',...
	'ToolTipString','clears lower results subpane of Graphics window',...
	'FontSize',Fs,'ForegroundColor','b',...
	'Callback',[...
		'spm_results_ui(''ClearPane''); ',...
		'spm_input(''!DeleteInputObj'')'],...
	'Interruptible','on','Enable','on',...
	'DeleteFcn','clear,clc,spm_clf(''Graphics'')',...
	'Position',[325 028 045 018].*WS);
hExit  = uicontrol(Finter,'Style','PushButton','String','exit',...
	'ToolTipString','exit the results section',...
	'FontSize',Fs,'ForegroundColor','r',...
	'Callback','spm_clf(''Interactive''), spm_clf(''Graphics''), clear',...
	'Interruptible','on','Enable','on',...
	'Position',[325 010 045 018].*WS);
hHelp  = uicontrol(Finter,'Style','PushButton','String','?',...
	'ToolTipString','results section help',...
	'FontSize',Fs,'ForegroundColor','g',...
	'Callback','spm_help(''spm_results.m'')',...
	'Interruptible','on','Enable','on',...
	'Position',[370 010 020 036].*WS);


%-Initialise registry in hReg frame object
%-----------------------------------------------------------------------
[hReg,xyz] = spm_XYZreg('InitReg',hReg,M,D,[0;0;0]);

%-Setup editable XYZ widgets & cross register with registry
%-----------------------------------------------------------------------
hFxyz      = spm_results_ui('DrawXYZgui',M,D,xyz,Finter);
spm_XYZreg('XReg',hReg,hFxyz,'spm_results_ui');

%-Set up buttons for Volume, Cluster & Voxel level utility functions
%-----------------------------------------------------------------------
[hFvolB,hFclusB,hFvoxB] = spm_results_ui('DrawButts',hReg,D,Finter,WS,Fs);

varargout  = {hReg};



%=======================================================================
case 'drawbutts'
%=======================================================================
% [hFvolB,hFclusB,hFvoxB] = spm_results_ui('DrawButts',hReg,D,Finter,WS,Fs)
if nargin<3, error('Insufficient arguments'), end
hReg = varargin{2};
D    = varargin{3};
if nargin<4,  Finter=spm_figure('FindWin','Interactive');
	else, Finter=varargin{4}; end
if nargin < 5, WS = spm('GetWinScale'); else,    WS = varargin{5}; end
if nargin < 6, Fs = 2*round(10*min(WS)/2); else, Fs = varargin{6}; end

if D(3) > 1, s3D = 'on'; else, s3D = 'off'; end


%-Volume level functions
%-----------------------------------------------------------------------
hFvolB = uicontrol(Finter,'Style','Frame','Position',[010 130 380 030].*WS);
uicontrol(Finter,'Style','Text','String','Volume level functions',...
	'Position',[020 152 125 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',Fs,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','p values','FontSize',Fs,...
	'ToolTipString','tabulate summary of p-values & statistics',...
	'Callback','spm_list(SPM,VOL)',...
	'Interruptible','on','Enable','on',...
	'Position',[015 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','multiplane','FontSize',Fs,...
	'ToolTipString','Launch MultiPlanar viewing toolbox',...
	'Tag','hB_MultiPlanar',...
	'Callback',...
		'hMP = spm_results_ui(''LaunchMP'',VOL.M,VOL.DIM,hReg,gcbo);',...
	'Interruptible','on','Enable','off',...
	'Position',[090 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','write','FontSize',Fs,...
	'ToolTipString','Write filtered SPM',...
	'Callback','spm_write_filtered(SPM,VOL)',...
	'Interruptible','on','Enable','on',...
	'Position',[165 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 135 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[315 135 070 020].*WS)

%-Cluster level functions
%-----------------------------------------------------------------------
hFclusB = uicontrol(Finter,'Style','Frame','Position',[010 090 380 030].*WS);
uicontrol(Finter,'Style','Text','String','Cluster level functions',...
	'Position',[020 112 125 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',Fs,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','maxima','FontSize',Fs,...
	'ToolTipString',['tabulate p-values & statistics for local maxima ',...
		'of current cluster'],...
	'Callback','spm_maxima(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[015 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','S.V.C.','FontSize',Fs,...
	'ToolTipString',['Small Volume Correction - corrected p-values ',...
		'for a small search region centered on this voxel'],...
	'Callback','spm_VOI(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[090 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[165 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 095 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
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
	'FontSize',Fs,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')
uicontrol(Finter,'Style','PushButton','String','plot','FontSize',Fs,...
	'ToolTipString','plotting of parameters & data at this voxel',...
	'Callback','spm_graph(SPM,VOL,DES,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[015 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','slices','FontSize',Fs,...
	'ToolTipString','overlay SPM on three slices of another image',...
	'Callback','spm_transverse(SPM,VOL,hReg)',...
	'Interruptible','on','Enable','on',...
	'Position',[090 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','sections','FontSize',Fs,...
	'ToolTipString',['overlay filtered SPM on orthogonal sections of',...
		'another image'],...
	'Callback','spm_sections(SPM,VOL,hReg)',...
	'Interruptible','on','Enable',s3D,...
	'Position',[165 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[240 055 070 020].*WS)
uicontrol(Finter,'Style','PushButton','String','','FontSize',Fs,...
	'ToolTipString','',...
	'Callback','',...
	'Interruptible','on','Enable','off',...
	'Position',[315 055 070 020].*WS)

varargout = {hFvolB,hFclusB,hFvoxB};

%=======================================================================
case 'drawxyzgui'
%=======================================================================
% hFxyz = spm_results_ui('DrawXYZgui',M,D,xyz,Finter)
if nargin<5,  Finter=spm_figure('FindWin','Interactive');
	else, Finter=varargin{5}; end
if nargin<4, xyz=[0;0;0]; else, xyz=varargin{4}; end
if nargin<3, error('Insufficient arguments'), end
D       = varargin{3};
M       = varargin{2};

xyz     = spm_XYZreg('RoundCoords',xyz,M,D);

%-Locate windows etc...
%-----------------------------------------------------------------------
WS      = spm('GetWinScale');
Fs      = 2*round(10*min(WS)/2);

%-Create & link control objects
%-----------------------------------------------------------------------
hFxyz = uicontrol(Finter,'Style','Frame','Position',[010 010 310 030].*WS);
uicontrol(Finter,'Style','Text','String','Co-ordinates',...
	'Position',[020 032 080 016].*WS,...
	'FontName','Times','FontWeight','Normal',...
	'FontSize',Fs,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w')


uicontrol(Finter,'Style','Text','String','x =',...
	'Position',[020 015 020 018].*WS,...
	'FontName','Times','FontSize',Fs,...
	'HorizontalAlignment','Center');
hX   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(1)),...
	'ToolTipString','enter x-coordinate',...
	'Position',[040 015 070 020].*WS,...
	'FontSize',Fs,'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hX',...
	'Callback','spm_results_ui(''EdWidCB'')');

uicontrol(Finter,'Style','Text','String','y =',...
	'Position',[120 015 020 018].*WS,...
	'FontName','Times','FontSize',Fs,...
	'HorizontalAlignment','Center')
hY   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(2)),...
	'ToolTipString','enter y-coordinate',...
	'Position',[140 015 070 020].*WS,...
	'FontSize',Fs,'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hY',...
	'Callback','spm_results_ui(''EdWidCB'')');

uicontrol(Finter,'Style','Text','String','z =',...
	'Position',[220 015 020 018].*WS,...
	'FontName','Times','FontSize',Fs,...
	'HorizontalAlignment','Center')
hZ   = uicontrol(Finter,'Style','Edit','String',sprintf('%.3f',xyz(3)),...
	'ToolTipString','enter z-coordinate',...
	'Position',[240 015 070 020].*WS,...
	'FontSize',Fs,'BackGroundColor',[.8,.8,1],...
	'HorizontalAlignment','Right',...
	'Tag','hZ',...
	'Callback','spm_results_ui(''EdWidCB'')');

set(hFxyz,'Tag','hFxyz','UserData',struct(...
	'hReg',	[],...
	'M',	M,...
	'D',	D,...
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
case 'edwidcb'		% Callback for editable widgets
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
	nxyz = spm_XYZreg('RoundCoords',nxyz,UD.M,UD.D);
end

if abs(xyz(d)-nxyz(d))>0
	UD.xyz = nxyz; set(hF,'UserData',UD)
	if ~isempty(UD.hReg), spm_XYZreg('SetCoords',nxyz,UD.hReg,hF); end
end

set(hC,'String',sprintf('%.3f',nxyz(d)))



%=======================================================================
case 'getcoords'	% Get current co-ordinates from XYZ widget
%=======================================================================
% xyz = spm_results_ui('GetCoords',hFxyz)
if nargin<2, hFxyz='Interactive'; else, hFxyz=varargin{2}; end
hFxyz     = spm_results_ui('FindXYZframe',hFxyz);
varargout = {getfield(get(hFxyz,'UserData'),'xyz')};



%=======================================================================
case 'setcoords'	% Set co-ordinates to XYZ widget
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
	[xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.D);
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
case 'findxyzframe'	% Find hFxyz frame
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
case 'plotui'
%=======================================================================
% spm_results_ui('PlotUi',hAx)
if nargin<2, hAx=gca; else, hAx=varargin{2}; end

WS = spm('GetWinScale');
Fs = 2*round(10*min(WS)/2);
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
	'FontSize',Fs,...
	'HorizontalAlignment','Left',...
	'ForegroundColor','w');

%-Controls
%-----------------------------------------------------------------------
h1 = uicontrol(Finter,'Style','CheckBox','String','Hold','FontSize',Fs,...
	'Value',strcmp(get(hAx,'NextPlot'),'add'),...
	'Callback',[...
		'if get(gcbo,''Value''), ',...
		    'set(get(gcbo,''UserData''),''NextPlot'',''add''), ',...
		'else, ',...
		    'set(get(gcbo,''UserData''),''NextPlot'',''replace''), ',...
		'end'],...
	'Interruptible','on','Enable','on',...
	'Position',[015 195 070 020].*WS);
h2 = uicontrol(Finter,'Style','CheckBox','String','Grid','FontSize',Fs,...
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
h3 = uicontrol(Finter,'Style','CheckBox','String','Box','FontSize',Fs,...
	'Value',strcmp(get(hAx,'Box'),'on'),...
	'Callback',[...
		'if get(gcbo,''Value''), ',...
		    'set(get(gcbo,''UserData''),''Box'',''on''), ',...
		'else, ',...
		    'set(get(gcbo,''UserData''),''Box'',''off''), ',...
		'end'],...
	'Interruptible','on','Enable','on',...
	'Position',[165 195 070 020].*WS);
h4 = uicontrol(Finter,'Style','PushButton','String','handle','FontSize',Fs,...
	'Callback','h=get(gcbo,''UserData'')',...
	'Interruptible','on','Enable','on',...
	'Position',[240 195 070 020].*WS);
h5 = uicontrol(Finter,'Style','PopUp','FontSize',Fs,...
	'String','Utils|Title|Xlabel|Ylabel|LineWidth|XLim|YLim',...
	'Callback','spm_results_ui(''PlotUiCB'')',...
	'Interruptible','off','Enable','on',...
	'Position',[315 195 070 020].*WS);

%-Handle storage for linking, and DeleteFcns for linked deletion
%-----------------------------------------------------------------------
set(hGraphUI,'UserData',[h1,h2,h3,h4,h5])
set([h1,h2,h3,h4,h5],'UserData',hAx)

set(hGraphUIbg,'UserData',...
	[hGraphUI,hGraphUIButtsF,hText,h1,h2,h3,h4,h5],...
	'DeleteFcn','delete(get(gcbo,''UserData''))')
set(hAx,'UserData',hGraphUIbg,...
	'DeleteFcn','delete(get(gcbo,''UserData''))')




%=======================================================================
case 'plotuicb'
%=======================================================================
% spm_results_ui('PlotUiCB')
hPM = gcbo;
v   = get(hPM,'Value');
set(hPM,'Value',1)

hAx = get(hPM,'UserData');
switch v
case 1	%-PullDown label
	return
case 2	%-Title
	h = get(hAx,'Title');
	set(h,'String',spm_input('Enter title:',-1,'s',get(h,'String')))
case 3	%-Xlabel
	h = get(hAx,'Xlabel');
	set(h,'String',spm_input('Enter X axis label:',-1,'s',get(h,'String')))
case 4	%-Ylabel
	h = get(hAx,'Ylabel');
	set(h,'String',spm_input('Enter Y axis label:',-1,'s',get(h,'String')))
case 5	%-LineWidth
	lw = spm_input('Enter LineWidth',-1,'e',get(hAx,'LineWidth'));
	set(hAx,'LineWidth',lw)
case 6	%-XLim
	XLim = spm_input('Enter XLim',-1,'e',get(hAx,'XLim'));
	set(hAx,'XLim',XLim)
case 7	%-YLim
	YLim = spm_input('Enter YLim',-1,'e',get(hAx,'YLim'));
	set(hAx,'YLim',YLim)
otherwise
	return
end



%=======================================================================
case 'clearpane'			%-Clear results subpane(s)
%=======================================================================
% Fgraph = spm_results_ui('ClearPane',F,RNP)
if nargin<3, RNP=0; else, RNP=1; end
if nargin<2, F=spm_figure('FindWin','Graphics'); else, F=varargin{2}; end
spm_input('!DeleteInputObj')
spm_figure('DeletePageControls',F)
set(0,'CurrentFigure',F)

kids = get(F,'Children');
for k=kids',
	if ishandle(k) & strcmp(get(k,'Type'),'axes'),
		un = get(k,'Units');
		set(k,'Units','normalized');
		pos = get(k,'Position');
		if pos(2) < 0.5,
			delete(k);
		else,
			set(k,'Units',un);
		end;
	end;
end;

% The original code:-
%subplot(2,1,2);
%if ~RNP | ~strcmp(get(gca,'NextPlot'),'add');
%	delete(gca); subplot(2,1,2); axis off;
%end



%=======================================================================
case 'launchmp'			%-Launch multiplanar toolbox
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
case 'delete'
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

