function varargout=spm(varargin)
% SPM: Statistical Parametric Mapping (startup function)
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  University College London
%
%_______________________________________________________________________
%
% SPM (Statistical Parametric Mapping) is a package for the analysis
% functional brain mapping experiments. It is the in-house package of
% the Wellcome Department of Cognitive Neurology, and is freely
% available to the scientific community.
% 
% Theoretical, computational and other details of the package are
% available in SPM's "Help" facility. This can be launched from the
% main SPM Menu window using the "Help" button, or directly from the
% command line using the command `spm_help`.
%
% Details of this release are availiable via the "About SPM" help topic
% (file spm.man), accessible from the SPM splash screen.  (Or type
% `spm_help spm.man` in the MatLab command window)
% 
% This spm function initialises the default parameters, and displays a
% splash screen with buttons leading to the PET(SPECT) & fMRI
% modalities Alternatively, `spm('pet')` and `spm('fmri')`
% (equivalently `spm pet` and `spm mri`) lead directly to the respective
% modality interfaces.
%
% Once the modality is chosen, (and it can be toggled mid-session) the
% SPM user interface is displayed. This provides a constant visual
% environment in which data analysis is implemented. The layout has
% been designed to be simple and at the same time show all the
% facilities that are available. The interface consists of three
% windows: A menu window with pushbuttons for the SPM routines (each
% button has a 'CallBack' string which launches the appropriate
% function/script); A blank panel used for interaction with the user;
% And a graphics figure with various editing and print facilities (see
% spm_figure.m). (These windows are 'Tag'ged 'Menu', 'Interactive', and
% 'Graphics' respectively, and should be referred to by their tags
% rather than their figure numbers.)
%
% Further interaction with the user is (mainly) via questioning in the
% 'Interactive' window (managed by spm_input), and file selection
% (managed by spm_get). See the help on spm_input.m and spm_get.m for
% details on using these functions.
%
% If a "message of the day" file named spm_motd.man exists in the SPM
% directory (alongside spm.m) then it is displayed in the Graphics
% window on startup.
%
% Arguments to this routine (spm.m) lead to various setup facilities,
% mainly of use to SPM power users and programmers. See programmers
% FORMAT & help in the main body of spm.m
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `spm Welcome` is equivalent to   )
%( `spm('Welcome')`.                                                   )
%
% FORMAT spm
% Defaults to spm('Welcome')
%
% FORMAT spm('Welcome')
% Clears command window, deletes all figures, prints welcome banner and
% splash screen, sets window defaults.
%
% FORMAT spm('AsciiWelcome')
% Prints ASCII welcome banner in MatLab command window.
%
% FORMAT spm('PET') spm('FMRI')
% Closes all windows and draws new Menu, Interactive, and Graphics
% windows for an SPM session. The buttons in the Menu window launch the
% main analysis routines.
%
% FORMAT Fmenu = spm('CreateMenuWin',Vis)
% Creates SPM menu window, 'Tag'ged 'Menu'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% Finter = FORMAT spm('CreateIntWin',Vis)
% Creates an SPM Interactive window, 'Tag'ged 'Interactive'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% FORMAT spm('ChMod',Modality)
% Changes modality of SPM: Currently SPM supports PET & MRI modalities,
% each of which have a slightly different Menu window and different
% defaults. This function switches to the specified modality, setting
% defaults and displaying the relevant buttons.
%
% FORMAT spm('defaults',Modality)
% Sets default global variables for the specified modality.
%
% FORMAT [Modality,ModNum]=spm('CheckModality',Modality)
% Checks the specified modality against those supported, returns
% upper(Modality) and the Modality number, it's position in the list of
% supported Modalities.
%
% FORMAT A=spm('GetWinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. A=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
%
% Rect = spm('WinSize',Win,raw)
% Returns sizes and positions for SPM windows.
% Win  - 'Menu', 'Interactive', 'Graphics', or '0'
%      -  Window whose position is required. Only first character is
%         examined. '0' returns size of root workspace.
% raw  - If specified, then positions are for a 1152 x 900 Sun display.
%        Otherwise the positions are scaled for the current display.
%
% FORMAT SPMdir=spm('Dir',Fname)
% Returns the directory containing the version of spm in use,
% identified as the first in MATLABPATH containing the Mfile spm (this
% file) (or Fname if specified).
%
% FORMAT SPMver=spm('Ver',Fname,ReDo)
% Returns the current version of SPM, identified by the top line of the
% Contents.m file in the directory containing the currently used file
% Fname (defaults on omission or empty to 'spm'). Inside SPM, the
% version is saved as the UserData of the Menu window, for quick
% retrieval. Specifying the ReDo argument forces re-evaluation.
%
% FORMAT [c,cName]=spm('Colour')
% Returns the RGB triple and a description for the current en-vogue SPM
% colour, the background colour for the Menu and Help windows.
%
% CmdLine = spm('isGCmdLine')
% Returns true if global CMDLINE exists and is itself true.
%
% FORMAT v = spm('MLver')
% Returns MatLab version, truncated to major & minor revision numbers
%
% FORMAT spm('SetCmdWinLabel',WinStripe,IconLabel)
% Sets the names on the headers and icons of Sun command tools.
% WinStripe defaults to a summary line identifying the user, host and
% MatLab version; IconLabel to 'MatLab'.
%
% FORMAT spm('UtilPullDownCB',h)
% Callback handler for "Utilities" PullDown menu in the SPM Menu window
%
% FORMAT User = spm('GetUser')
% Returns current user, culled from the USER environment variable
%
%_______________________________________________________________________

%-Parameters
%-----------------------------------------------------------------------
Modalities = str2mat('PET','FMRI');

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Welcome'; else, Action = varargin{1}; end


switch lower(Action), case 'welcome'
%=======================================================================

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
S = get(0,'ScreenSize');
F = figure('IntegerHandle','off',...
	'Name',[spm('GetUser'),' - SPM'],'NumberTitle','off',...
	'Tag','Welcome',...
	'Position',[S(3)/2-300,S(4)/2-140,500,280],...
	'Resize','off',...
	'Pointer','Watch',...
	'Color',[1 1 1]*.8,...
	'MenuBar','none',...
	'HandleVisibility','off',...
	'Visible','off');

%-Text
%-----------------------------------------------------------------------
hA = axes('Parent',F,'Position',[0 0 100/500 280/280],'Visible','Off');
text(0.5,0.5,'SPM',...
	'Parent',hA,...
	'FontName','Times','FontSize',96,...
	'FontAngle','Italic','FontWeight','Bold',...
	'Rotation',90,...
	'VerticalAlignment','Middle','HorizontalAlignment','Center',...
	'Color',[1 1 1]*.6);

% uicontrol(F,'Style','Text','Position',[000 190 110 100],...
% 	'String','S','FontSize',84,...
% 	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
% 	'HorizontalAlignment','center',...
% 	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);
% uicontrol(F,'Style','Text','Position',[000 095 110 100],...
% 	'String','P','FontSize',84,...
% 	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
% 	'HorizontalAlignment','center',...
% 	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);
% uicontrol(F,'Style','Text','Position',[000 000 110 100],...
% 	'String','M','FontSize',84,...
% 	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
% 	'HorizontalAlignment','center',...
% 	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);


uicontrol(F,'Style','Text','Position',[080 245 420 030],...
	'String','Statistical Parametric Mapping',...
	'FontName','Times','FontSize',18,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);


uicontrol(F,'Style','Frame','Position',[110 120 380 120]);
uicontrol(F,'Style','Frame','Position',[110 020 380 090]);

uicontrol(F,'Style','Text','String',spm('Ver'),...
	'Position',[112 200 376 030],...
	'FontName','Times','FontSize',14,'FontWeight','Bold',...
	'ForegroundColor','b')

uicontrol(F,'Style','Text','Position',[112 180 376 020],'String',...
	'The Wellcome Department of Cognitive Neurology',...
	'FontName','Times','FontSize',12,'FontWeight','Bold')
uicontrol(F,'Style','Text', 'Position',[112 160 376 020],...
	'FontName','Times','FontSize',12,'String',...
	'The Institute of Neurology')
uicontrol(F,'Style','Text', 'Position',[112 140 376 020],...
	'FontName','Times','FontSize',12,'String',...
	'University College London')

%-Objects with Callbacks - PET, fMRI, About SPM, SPMweb
%-----------------------------------------------------------------------
set(F,'DefaultUicontrolFontSize',12,'DefaultUicontrolInterruptible','on')
uicontrol(F,'String','PET and SPECT',...
	'Position',[140 066 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''PET'')',...
	'ForegroundColor',[0 1 1])
uicontrol(F,'String','fMRI time-series',...
	'Position',[310 066 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''FMRI'')',...
	'ForegroundColor',[0 1 1])
uicontrol(F,'String','About SPM',...
	'Position',[140 030 100 030],...
	'CallBack','spm_help(''spm.man'')',...
	'ForegroundColor','g')
uicontrol(F,'String','SPMweb',...
	'FontWeight','Bold','FontName','Courier',...
	'Position',[250 030 100 030],...
	'CallBack',['set(gcbf,''Pointer'',''Watch''),',...
			'web(''http://www.fil.ion.ucl.ac.uk/spm'');',...
			'set(gcbf,''Pointer'',''Arrow'')'],...
	'ForegroundColor','k')
uicontrol(F,'String','Quit',...
	'Position',[360 030 100 030],...
	'CallBack','delete(gcbf)',...
	'Interruptible','off',...
	'ForegroundColor','r')

set(F,'Pointer','Arrow','Visible','on')


case 'asciiwelcome'
%=======================================================================
disp( ' ___  ____  __  __                                                  ')
disp( '/ __)(  _ \(  \/  )  Statistical Parametric Mapping                 ')
disp( '\__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology ')
disp(['(___/(__)  (_/\/\_)  ',spm('Ver'),' - http://www.fil.ion.ucl.ac.uk/spm'])
fprintf('\n')


case {'pet','fmri'}
%=======================================================================
% spm(Modality)

%-Turn on warning messages for ML5 debugging
warning always, warning backtrace

%-Initialisation and workspace canonicalisation
%-----------------------------------------------------------------------
clc, spm('SetCmdWinLabel')
spm('AsciiWelcome'),			fprintf('\n\nInitialising SPM')
Modality = upper(Action);					fprintf('.')
delete(get(0,'Children')),					fprintf('.')

%-Draw SPM windows
%-----------------------------------------------------------------------
Fmenu = spm('CreateMenuWin','off');				fprintf('.')
Finter = spm('CreateIntWin','off');				fprintf('.')
spm_figure('WaterMark',Finter,spm('Ver'),'',45),		fprintf('.')
Fgraph = spm_figure('Create','Graphics','Graphics','off');	fprintf('.')

Fmotd = [spm('Dir'),'/spm_motd.man'];
if exist(Fmotd), spm_help('!Disp',Fmotd,'',Fgraph,spm('Ver')); end
								fprintf('.')

%-Load startup global defaults
%-----------------------------------------------------------------------
spm_defaults,							fprintf('.')

%-Setup for current modality
%-----------------------------------------------------------------------
spm('ChMod',Modality),						fprintf('.')

%-Reveal windows
%-----------------------------------------------------------------------
set([Fmenu,Finter,Fgraph],'Visible','on')
							fprintf('done\n\n')

%-Print present working directory
%-----------------------------------------------------------------------
fprintf('SPM present working directory:\n\t%s\n',pwd)




case 'createmenuwin'
%=======================================================================
% Fmenu = spm('CreateMenuWin',Vis)
if nargin<2, Vis='on'; else, Vis=varargin{2}; end

%-Close any existing 'Menu' 'Tag'ged windows
close(spm_figure('FindWin','Menu'))

%-Get size and scalings and create Menu window
%-----------------------------------------------------------------------
WS   = spm('GetWinScale');
Rect = spm('WinSize','Menu','raw').*WS;

SPMver = spm('Ver','',1);

Fmenu = figure('IntegerHandle','off',...
	'Name',[spm('GetUser'),' - ',SPMver],'NumberTitle','off',...
	'Tag','Menu',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'UserData',SPMver,...
	'MenuBar','none',...
	'DefaultUicontrolFontSize',2*round(12*min(WS)/2),...
	'DefaultUicontrolInterruptible','on',...
	'Visible','off');

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 145 380 295].*WS,'Tag','Empty')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 320 360 110].*WS,'Tag','Empty')
uicontrol(Fmenu,'Style','Text','String','Spatial',...
	'Position',[100 400 200 20].*WS,'Tag','Empty',...
	'ForegroundColor','w')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 235 360 075].*WS,'Tag','Empty');
	uicontrol(Fmenu,'Style','Text','String','Analysis',...
	'Position',[100 285 200 20].*WS,'Tag','Empty',...
	'ForegroundColor','w')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 155 360 070].*WS,'Tag','Empty');
uicontrol(Fmenu,'Style','Text','String','Results',...
	'Position',[100 200 200 20].*WS,'Tag','Empty',...
	'ForegroundColor','w' )

uicontrol(Fmenu,'Style','Text',...
	'String','SPM for PET/SPECT',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','PET','Visible','off')
uicontrol(Fmenu,'Style','Text',...
	'String','SPM for functional MRI',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','FMRI','Visible','off')

uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 010 380 112].*WS,'Tag','Empty');

%-Objects with Callbacks - main spm_*_ui.m routines
%-----------------------------------------------------------------------

%-Spatial
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Realign',	'Position',[040 370 080 30].*WS,...
	'CallBack','spm_realign')

uicontrol(Fmenu,'String','Normalize',	'Position',[150 350 100 30].*WS,...
	'CallBack','spm_sn3d')

uicontrol(Fmenu,'String','Smooth',	'Position',[280 370 080 30].*WS,...
	'CallBack','spm_smooth_ui')

uicontrol(Fmenu,'String','Coregister',	'Position',[040 330 080 30].*WS,...
	'CallBack','spm_coregister')

uicontrol(Fmenu,'String','Segment',	'Position',[280 330 080 30].*WS,...
	'CallBack','spm_segment')

%-Statistics
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Statistics',	'Position',[040 245 140 30].*WS,...
	'CallBack','spm_spm_ui',...
	'Visible','off',		'Tag','PET')

uicontrol(Fmenu,'String','Statistics',	'Position',[040 245 140 30].*WS,...
	'CallBack','spm_fmri_spm_ui',...
	'Visible','off',		'Tag','FMRI')

uicontrol(Fmenu,'String','Eigenimages',	'Position',[220 245 140 30].*WS,...
	'CallBack','spm_svd_ui')

%-Results
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','SPM{F}',	'Position',[045 165 070 30].*WS,...
	'CallBack','spm_projectionsF_ui;')

uicontrol(Fmenu,'String','Results',	'Position',[165 165 070 30].*WS,...
	'CallBack','spm_results')

uicontrol(Fmenu,'String','SPM{Z}',	'Position',[285 165 070 30].*WS,...
	'CallBack','spm_projections_ui;')

%-Utility buttons (first line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Display',	'Position',[020 088 082 024].*WS,...
	'CallBack','spm_image')

uicontrol(Fmenu,'String','',		'Position',[112 088 083 024].*WS,...
	'CallBack','',			'Enable','off')

uicontrol(Fmenu,'String','Render',	'Position',[205 088 083 024].*WS,...
	'CallBack','spm_render')

uicontrol(Fmenu,'Style','PopUp','String',Modalities,...
	'Tag','Modality',		'Position',[298 088 082 024].*WS,...
	'CallBack',...
	[	'if get(gco,''Value'')~=get(gco,''UserData''),',...
			'spm(''ChMod'',get(gco,''Value'')),',...
		'end'],...
					'Interruptible','off')

%-Utility buttons (second line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','CD',		'Position',[020 054 082 024].*WS,...
	'CallBack',[...
		'global CWD,',...
		'CWD=spm_get(-1,[],''Select new working directory'');',...
		'cd(CWD), clear CWD,',...
		'fprintf(''\nSPM working directory:\n\t%s\n\n>> '',pwd)'])

uicontrol(Fmenu,'String','Mean',	'Position',[112 054 083 024].*WS,...
	'CallBack','spm_average')

uicontrol(Fmenu,'String','ImCalc',	'Position',[205 054 083 024].*WS,...
	'CallBack','spm_image_funks')

uicontrol(Fmenu,'String','HDR edit',	'Position',[298 054 082 024].*WS,...
	'CallBack','spm_header_edit')

%-Utility buttons (third line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Help',	'Position',[020 020 082 024].*WS,...
	'CallBack','spm_help',...
	'ForeGroundColor','g')

uicontrol(Fmenu,'Style','PopUp',...
	'String','Utilities|Analyze|GhostView|Run mFile',...
	'Tag','UtilPullDown',		'Position',[112 020 083 024].*WS,...
	'CallBack','spm(''UtilPullDownCB'',gcbo)',...
	'UserData',{	'!analyze',...
			['unix([''ghostview '',spm_get(1,''.ps'',''Select ',...
				'PostScript file to view''),'' &'']);'],...
			'run(spm_get(1,''*.m'',''Select mFile to run''))'} )

uicontrol(Fmenu,'String','Defaults',	'Position',[205 020 083 024].*WS,...
	'CallBack','spm_defaults_edit')

uicontrol(Fmenu,'String','Quit',	'Position',[298 020 082 024].*WS,...
	'ForeGroundColor','r',		'Interruptible','off',...
	'CallBack',[...
		'delete(get(0,''Children'')),',...
		'clear all,clc,fprintf(''Bye...\n\n>> '')'])
%-----------------------------------------------------------------------
set(Fmenu,'CloseRequestFcn',[...
		'delete(get(0,''Children'')),',...
		'clear all,clc,fprintf(''Bye...\n\n>> '')'])

set(Fmenu,'Visible',Vis)
varargout = {Fmenu};



case 'createintwin'
%=======================================================================
% Finter = spm('CreateIntWin',Vis)
if nargin<2, Vis='on'; else, Vis=varargin{2}; end

%-Close any existing 'Interactive' 'Tag'ged windows
close(spm_figure('FindWin','Interactive'))

%-Create SPM Interactive window
Rect = spm('WinSize','Interactive');
Finter = figure('IntegerHandle','off',...
	'Tag','Interactive',...
	'Name','','NumberTitle','off',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.7,...
	'MenuBar','none',...
	'DefaultUicontrolFontSize',2*round(12*min(spm('GetWinScale'))/2),...
	'DefaultUicontrolInterruptible','on',...
	'Visible',Vis);
varargout = {Finter};


case 'chmod'
%=======================================================================
% spm('ChMod',Modality)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality = ''; else, Modality = varargin{2}; end
[Modality,ModNum] = spm('CheckModality',Modality);

if strcmp(Modality,'PET'), OModality = 'FMRI'; else, OModality='PET'; end


%-Sort out global defaults
%-----------------------------------------------------------------------
spm('defaults',Modality)

%-Sort out visability of appropriate controls on Menu window
%-----------------------------------------------------------------------
Fmenu = spm_figure('FindWin','Menu');
if isempty(Fmenu), error('SPM Menu window not found'), end

set(findobj(Fmenu,'Tag',OModality),'Visible','off')
set(findobj(Fmenu,'Tag', Modality),'Visible','on')
set(findobj(Fmenu,'Tag','Modality'),'Value',ModNum,'UserData',ModNum)

%-For fMRI, check descriptors limit
%-----------------------------------------------------------------------
if strcmp(Modality,'FMRI')
	[s,w] = unix('limit');
	d     = findstr(w,'descriptors');
	w     = eval(w((d + 11):length(w)));
	if (w < 256) warndlg(['To increase; quit SPM & MatLab, type ',...
		'''unlimit'', and restart'],...
		sprintf('WARNING: file descriptors = %d',w))
	end
end


case 'defaults'
%=======================================================================
% spm('defaults',Modality)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=varargin{2}; end
Modality = spm('CheckModality',Modality);

%-Set MODALITY
%-----------------------------------------------------------------------
% clear global
global MODALITY
MODALITY = Modality;

%-Set global defaults (global variables)
%-----------------------------------------------------------------------
global SWD TWD DESCRIP
global CWD

SWD	 = spm('Dir');					% SPM directory
CWD	 = pwd;						% Working directory
TWD	 = getenv('SPMTMP');				% Temporary directory
if isempty(TWD)
	TWD = '/tmp';
end

%-Get global modality defaults
%-----------------------------------------------------------------------
global PET_UFp PET_DIM PET_VOX PET_TYPE PET_SCALE PET_OFFSET PET_ORIGIN PET_DESCRIP
global fMRI_UFp fMRI_DIM fMRI_VOX fMRI_TYPE fMRI_SCALE fMRI_OFFSET fMRI_ORIGIN fMRI_DESCRIP

%-Load startup defaults (if not already done so)
%-----------------------------------------------------------------------
if isempty(PET_DIM), spm_defaults, end

%-Set Modality specific default (global variables)
%-----------------------------------------------------------------------
global UFp DIM VOX SCALE TYPE OFFSET ORIGIN
if strcmp(Modality,'PET')
	DIM	= PET_DIM;				% Dimensions [x y z]
	VOX	= PET_VOX;				% Voxel size [x y z]
	SCALE	= PET_SCALE;				% Scaling coeficient
	TYPE	= PET_TYPE;				% Data type
	OFFSET	= PET_OFFSET;		    		% Offset in bytes
	ORIGIN	= PET_ORIGIN;				% Origin in voxels
	UFp	= PET_UFp;				% Upper tail F-prob
elseif strcmp(Modality,'FMRI')
	DIM	= fMRI_DIM;				% Dimensions {x y z]
	VOX	= fMRI_VOX;				% Voxel size {x y z]
	SCALE	= fMRI_SCALE;				% Scaling coeficient
	TYPE	= fMRI_TYPE;				% Data type
	OFFSET	= fMRI_OFFSET;	   			% Offset in bytes
	ORIGIN	= fMRI_ORIGIN;				% Origin in voxels
	UFp	= fMRI_UFp;				% Upper tail F-prob

elseif strcmp(Modality,'UNKNOWN')
else
	error('Illegal Modality')
end


case 'checkmodality'
%=======================================================================
% [Modality,ModNum] = spm('CheckModality',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=upper(varargin{2}); end
if isempty(Modality)
	global MODALITY
	if ~isempty(MODALITY); Modality = MODALITY;
	else, Modality = 'UNKNOWN'; end
end
if isstr(Modality)
	ModNum = find(all(Modalities(:,1:length(Modality))'==...
			Modality'*ones(1,size(Modalities,1))));
else
	if ~any(Modality==[1:size(Modalities,1)])
		Modality = 'ERROR';
		ModNum   = [];
	else
		ModNum   = Modality;
		Modality = deblank(Modalities(ModNum,:));
	end
end

if isempty(ModNum), error('Unknown Modality'), end
varargout = {upper(Modality),ModNum};


case 'getwinscale'
%=======================================================================
% spm('GetWinScale')
S   = get(0,'ScreenSize');
varargout = {[S(3)/1152 S(4)/900 S(3)/1152 S(4)/900]};


case 'winsize'
%=======================================================================
% Rect = spm('WinSize',Win,raw)
if nargin<3, raw=0; else, raw=1; end
if nargin<2, Win=''; else, Win=varargin{2}; end

Rect = [	[108 429 400 445];...
		[108 008 400 395];...
		[515 008 600 865] ];

WS = spm('GetWinScale');

if isempty(Win)
	WS = ones(3,1)*WS;
elseif upper(Win(1))=='M'
	%-Menu window
	Rect = Rect(1,:);
elseif upper(Win(1))=='I'
	%-Interactive window
	Rect = Rect(2,:);
elseif upper(Win(1))=='G'
	%-Graphics window
	Rect = Rect(3,:);
elseif Win(1)=='0'
	%-Root workspace
	Rect = get(0,'ScreenSize');
else
	error('Unknown Win type');
end

if ~raw, Rect = Rect.*WS; end
varargout = {Rect};


case 'dir'
%=======================================================================
% spm('Dir',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='spm'; else, Mfile=varargin{2}; end
SPMdir = which(Mfile);
if ~isstr(SPMdir)
	error(['Can''t find ',Mfile,' on MATLABPATH']); end
tmp    = max(find(SPMdir=='/'))-1;
if tmp, SPMdir = SPMdir(1:tmp); end
varargout = {SPMdir};


case 'ver'
%=======================================================================
% spm('Ver',Mfile,ReDo)
if nargin<3, ReDo=0; else, ReDo=1; end
if nargin<2, Mfile=''; else, Mfile=varargin{2}; end
if isempty(Mfile), Mfile='spm'; end
Fmenu = spm_figure('FindWin','Menu');

%-See if SPM Menu window exists - It's UserData should contain the version
%-----------------------------------------------------------------------
if ~ReDo
	if ~isempty(Fmenu)
		SPMver = get(Fmenu,'UserData');
		if ~isempty(SPMver), varargout={SPMver}; return, end
	end
end

%-Work version out from file
%-----------------------------------------------------------------------
SPMdir = spm('Dir',Mfile);
CFile  = [SPMdir,'/Contents.m'];
if exist(CFile)
	fid  = fopen(CFile,'r');
	SPMver = setstr([fread(fid,80,'char')',setstr(10)]);
	fclose(fid);
	SPMver(1:max(1,min(find(SPMver~='%' & SPMver~=' '))-1))=[];
	SPMver = SPMver(1:min(find(SPMver==10))-1);
else
	SPMver = 'SPM';
end

%-Store version in UserData of SPM Menu window, if it exists
%-----------------------------------------------------------------------
if ~isempty(Fmenu), set(Fmenu,'UserData',SPMver); end

varargout = {SPMver};


case 'colour'
%=======================================================================
% spm('Colour')
%-----------------------------------------------------------------------
%-Developmental livery
varargout = {[0.7,1.0,0.7],'Lime Green'};
%-Distribution livery
% varargout = {[0.8 0.8 1.0],'Diluted Blackcurrent Purple'};


case 'isgcmdline'
%=======================================================================
% CmdLine = spm('isGCmdLine')
global CMDLINE
if isempty(CMDLINE), varargout = {0}; else, varargout = {CMDLINE}; end


case 'mlver'
%=======================================================================
% v = spm('MLver')
v = version; tmp = find(v=='.');
if length(tmp)>1, varargout={v(1:tmp(2)-1)}; end


case 'setcmdwinlabel'
%=======================================================================
% spm('SetCmdWinLabel',WinStripe,IconLabel)
%-----------------------------------------------------------------------

%-Only label Sun command tools
%-----------------------------------------------------------------------
Term        = getenv('TERM');
if ~strcmp(Term,'sun-cmd'), return, end

%-Work out label text
%-----------------------------------------------------------------------
User        = spm('GetUser');
[null,Host] = unix('echo `hostname` | sed -e ''s/\..*$//''');
Host        = Host(1:length(Host)-1);
v           = spm('MLver');

if nargin<3, IconLabel = ['MatLab',v(1)]; end
if nargin<2, WinStripe = [User,' - ',Host,' : MatLab ',v]; end

%-Set window stripe
%-----------------------------------------------------------------------
disp([']l' WinStripe '\]L' IconLabel '\'])


case 'utilpulldowncb'
%=======================================================================
% spm('UtilPullDownCB',h)
if nargin<2, h=gcbo; else, h=varargin{2}; end
v   = get(h,'Value');
if v==1, return, end
set(h,'Value',1)
CBs = get(h,'UserData');
evalin('base',CBs{v-1})


case 'getuser'
%=======================================================================
% User = spm('GetUser')
User = getenv('USER');
if isempty(User), User='user'; end
varargout = {User};

case 'beep'
%=======================================================================
% spm('Beep')
fprintf('%c',7)


otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end
