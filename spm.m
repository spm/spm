function varargout=spm(varargin)
% SPM: Statistical Parametric Mapping (startup function)
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  Institute of Neurology, University College London
%
%_______________________________________________________________________
%
% SPM (Statistical Parametric Mapping) is a package for the analysis
% functional brain mapping experiments. It is the in-house package of
% the Wellcome Department of Cognitive Neurology, and is available to
% the scientific community as freeware.
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
% FORMAT WS=spm('GetWinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. WS=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
% (Function duplicated in spm_figure.m, repeated to reduce inter-dependencies.)
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
% FORMAT SPMver=spm('Ver',Fname,ReDo,Cache)
% Returns the current version of SPM, identified by the top line of the
% Contents.m file in the directory containing the currently used file
% Fname (defaults on omission or empty to 'spm'). Inside SPM, the
% version can be saved in a structure as the UserData of the Menu
% window, for quick retrieval. If Cache [default 1] is true, then the
% version information is so cached. If ReDo [default 0] is true this
% forces re-evaluation, regardless of the presence of the information
% in the MenuWin's UserData. (For toolboxes, should therefore use
% spm('Ver',Fname,1,0) where Fname is one of the toolbox routines.)
%
% FORMAT [c,cName]=spm('Colour')
% Returns the RGB triple and a description for the current en-vogue SPM
% colour, the background colour for the Menu and Help windows.
%
% FORMAT [v1,v2,...] = spm('GetGlobal',name1,name2,...)
% Returns values of global variables (without declaring them global)
% name1, name2,... - name strings of desired globals
% a1, a2,...       - corresponding values of global variables with given names
%                    ([] is returned as value if global variable doesn't exist)
%
% FORMAT CmdLine = spm('isGCmdLine')
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
% FORMAT spm('PopUpCB',h)
% Callback handler for PopUp UI menus with multiple callbacks as cellstr UserData
%
% FORMAT User = spm('GetUser')
% Returns current user, culled from the USER environment variable
%
% FORMAT spm('Beep')
% plays the keyboard beep!
%
% FORMAT spm('time')
% Returns the current time and date as hh:mm dd/mm/yyyy
%
% FORMAT spm('Pointer',Pointer)
% Changes pointer on all SPM (HandleVisible) windows to type Pointer
% Pointer defaults to 'Arrow'. Robust to absence of windows
%
% FORMAT SPMid = spm('FnBanner',Fn,FnV)
% Prints a function start banner, for version FnV of function Fn, & datestamps
% Fn    - Function name (string)
% FnV   - Function version (string)
% SPMid - ID string: [SPMver: Fn (FnV)] 
%
% FORMAT [Finter,Fgraph,CmdLine] = spm('FnUIsetup',Iname,bGX,CmdLine)
% Robust UIsetup procedure for functions:
%   Returns handles of 'Interactive' and 'Graphics' figures.
%   Creates 'Interactive' figure if ~CmdLine, creates 'Graphics' figure if bGX.
% Iname   - Name for 'Interactive' window
% bGX     - Need a Graphics window? [default 1]
% CmdLine - CommandLine usage? [default spm('isGCmdLine')]
% Finter  - handle of 'Interactive' figure
% Fgraph  - handle of 'Graphics' figure
% CmdLine - CommandLine usage?
%
% FORMAT F = spm('FigName',Iname,F,CmdLine)
% Set name of figure F to [spm('GetUser'),' - ',Iname] if ~CmdLine
% Robust to absence of figure.
% Iname      - Name for figure
% F (input)  - Handle (or 'Tag') of figure to name [default 'Interactive']
% CmdLine    - CommandLine usage? [default spm('isGCmdLine')]
% F (output) - Handle of figure named
%
% FORMAT Fs = spm('Show')
% Opens all SPM figure windows (with HandleVisibility) using `figure`.
%   Maintains current figure.
% Fs - vector containing all HandleVisible figures (i.e. get(0,'Children'))
%
% FORMAT spm('Clear',Finter, Fgraph)
% Clears and resets SPM-GUI, clears and timestamps MatLab command window
% Finter  - handle or 'Tag' of 'Interactive' figure [default 'Interactive']
% Fgraph  - handle or 'Tag' of 'Graphics' figure [default 'Graphics']
%
% FORMAT spm('Help',varargin)
% Merely a gateway to spm_help(varargin) - so you can type "spm help"
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
if all(S==1), error('Can''t open any graphics windows...'), end
[SPMver,SPMc] = spm('Ver','',1);

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

uicontrol(F,'Style','Text','Position',[080 245 420 030],...
	'String','Statistical Parametric Mapping',...
	'FontName','Times','FontSize',18,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);


uicontrol(F,'Style','Frame','Position',[110 150 380 095]);
uicontrol(F,'Style','Frame','Position',[110 015 380 087]);

uicontrol(F,'Style','Text','String',SPMver,...
	'ToolTipString','by the FIL methods group',...
	'Position',[112 200 376 030],...
	'FontName','Times','FontSize',14,'FontWeight','Bold',...
	'ForegroundColor','b')

uicontrol(F,'Style','Text','Position',[112 180 376 020],...
	'String','The Wellcome Department of Cognitive Neurology',...
	'ToolTipString','',...
	'FontName','Times','FontSize',12,'FontWeight','Bold')
uicontrol(F,'Style','Text', 'Position',[112 160 376 020],...
	'String','The Institute of Neurology, University College London',...
	'ToolTipString','',...
	'FontName','Times','FontSize',12)

uicontrol(F,'Style','Text','String',SPMc,'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontSize',8,...
	'HorizontalAlignment','center',...
	'Position',[110 003 380 010])

%-Objects with Callbacks - PET, fMRI, About SPM, SPMweb
%-----------------------------------------------------------------------
set(F,'DefaultUicontrolFontSize',12,'DefaultUicontrolInterruptible','on')
uicontrol(F,'String','PET and SPECT',...
	'ToolTipString','launch SPM-GUI in PET/SPECT modality',...
	'Position',[140 061 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''PET'')',...
	'ForegroundColor',[0 1 1])
uicontrol(F,'String','fMRI time-series',...
	'ToolTipString','launch SPM-GUI in fMRI modality',...
	'Position',[310 061 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''FMRI'')',...
	'ForegroundColor',[0 1 1])
uicontrol(F,'String','About SPM',...
	'ToolTipString','launch SPMhelp browser - about SPM',...
	'Position',[140 025 100 030],...
	'CallBack','spm_help(''spm.man'')',...
	'ForegroundColor','g')
uicontrol(F,'String','SPMweb',...
	'FontWeight','Bold','FontName','Courier',...
	'ToolTipString',...
	'launch web browser - http://www.fil.ion.ucl.ac.uk/spm',...
	'Position',[250 025 100 030],...
	'CallBack',['set(gcbf,''Pointer'',''Watch''),',...
			'web(''http://www.fil.ion.ucl.ac.uk/spm'');',...
			'set(gcbf,''Pointer'',''Arrow'')'],...
	'ForegroundColor','k')
uicontrol(F,'String','Quit',...
	'ToolTipString','close this splash screen',...
	'Position',[360 025 100 030],...
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
delete(spm_figure('FindWin','Menu'))

%-Get size and scalings and create Menu window
%-----------------------------------------------------------------------
WS   = spm('GetWinScale');			%-Window scaling factors
Rect = spm('WinSize','Menu','raw').*WS;		%-Menu window rectangle
FS3  = 2*round(12*min(WS)/2);			%-Standard font size
FS2  = round(9*min(WS));			%-Smaller font size

[SPMver,SPMc] = spm('Ver','',1,1);

Fmenu = figure('IntegerHandle','off',...
	'Name',[spm('GetUser'),' - ',SPMver],'NumberTitle','off',...
	'Tag','Menu',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'UserData',struct('SPMver',SPMver,'SPMc',SPMc),...
	'MenuBar','none',...
	'DefaultUicontrolFontSize',FS3,...
	'DefaultUicontrolInterruptible','on',...
	'Visible','off');


%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 145 380 295].*WS)

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 320 360 110].*WS)
uicontrol(Fmenu,'Style','Text','String','spatial',...
	'Position',[025 405 350 020].*WS,...
	'ForegroundColor','w','FontName','Times','FontAngle','Italic')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 155 360 155].*WS)
uicontrol(Fmenu,'Style','Text','String','statistical',...
	'Position',[025 285 350 020].*WS,...
	'ForegroundColor','w','FontName','Times','FontAngle','Italic')

uicontrol(Fmenu,'Style','Text',...
	'String','SPM for PET/SPECT',...
	'ToolTipString','modality & defaults set for PET/SPECT',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','PET','Visible','off')
uicontrol(Fmenu,'Style','Text',...
	'String','SPM for functional MRI',...
	'ToolTipString','modality & defaults set for fMRI',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontAngle','Italic','FontWeight','Bold',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','FMRI','Visible','off')

uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 010 380 112].*WS);

uicontrol(Fmenu,'Style','Text','String',SPMc,'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName','Times','FontSize',6,...
	'HorizontalAlignment','center',...
	'Position',[020 002 360 008].*WS)

%-Objects with Callbacks - main spm_*_ui.m routines
%=======================================================================

%-Spatial
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Realign',	'Position',[040 370 080 030].*WS,...
	'ToolTipString','realignment',...
	'CallBack','spm_realign')

uicontrol(Fmenu,'String','Normalize',	'Position',[150 350 100 030].*WS,...
	'ToolTipString','spatial normalisation',...
	'CallBack','spm_sn3d')

uicontrol(Fmenu,'String','Smooth',	'Position',[280 370 080 030].*WS,...
	'ToolTipString','spatial smoothing with Gaussian kernel',...
	'CallBack','spm_smooth_ui')

uicontrol(Fmenu,'String','Coregister',	'Position',[040 330 080 030].*WS,...
	'ToolTipString','co-register images from disparate modalities',...
	'CallBack','spm_coregister')

uicontrol(Fmenu,'String','Segment',	'Position',[280 330 080 030].*WS,...
	'ToolTipString','segment',...
	'CallBack','spm_segment')

%-Statistical
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','PET/SPECT models','Position',[035 255 160 030].*WS,...
	'ToolTipString','general linear model setup for PET/SPECT',...
	'CallBack','spm_spm_ui(''cfg'',spm_spm_ui(''DesDefs_SPM98PET''))',...
	'Visible','off',		'Tag','PET'	,'Enable','on')
uicontrol(Fmenu,'String','fMRI models',	'Position',[035 255 160 030].*WS,...
	'ToolTipString',['general linear model setup & stats for serially ',...
		'correlated fMRI time series'],...
	'CallBack','spm_fmri_spm_ui',...
	'Visible','off',		'Tag','FMRI')
uicontrol(Fmenu,'String','Basic models','Position',[205 255 160 030].*WS,...
	'ToolTipString','basic stats  models for independent data',...
	'CallBack','spm_spm_ui(''cfg'',spm_spm_ui(''DesDefs_Stats''))')

uicontrol(Fmenu,'String','Results',	'Position',[035 215 330 030].*WS,...
	'ToolTipString','interrogate model, contrasts etc.',...
	'CallBack','spm_results')
uicontrol(Fmenu,'String','',		'Position',[035 215 160 030].*WS,...
	'ToolTipString','',...
	'CallBack','',	'Visible','off',	'Tag',''	,'Enable','off')
uicontrol(Fmenu,'String','',		'Position',[205 215 160 030].*WS,...
	'ToolTipString','',...
	'CallBack','',	'Visible','off',	'Tag',''	,'Enable','off')

uicontrol(Fmenu,'String','',		'Position',[035 170 160 030].*WS,...
	'ToolTipString','',...
	'CallBack','',	'Visible','on',		'Tag',''	,'Enable','off')

uicontrol(Fmenu,'String','Eigenimages',	'Position',[205 170 160 030].*WS,...
	'ToolTipString','PCA of adjusted data (not coded for SPM98 yet)',...
	'CallBack','spm_svd_ui'				,'Enable','off')


%-Utility buttons (first line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Display',	'Position',[020 088 082 024].*WS,...
	'ToolTipString','orthogonal sections',...
	'FontSize',FS2,			'CallBack','spm_image')

uicontrol(Fmenu,'String','Check Reg',	'Position',[112 088 083 024].*WS,...
	'ToolTipString','check image registration',...
	'FontSize',FS2,		'CallBack','spm_check_registration')

uicontrol(Fmenu,'String','Render',	'Position',[205 088 083 024].*WS,...
	'ToolTipString','render',...
	'FontSize',FS2,			'CallBack','spm_render')

uicontrol(Fmenu,'Style','PopUp','String',Modalities,...
	'ToolTipString','change modality PET<->fMRI',...
	'Tag','Modality',		'Position',[298 088 082 024].*WS,...
	'CallBack',[...
		'if isempty(get(gco,''UserData'')) | ',...
		'get(gco,''Value'')~=get(gco,''UserData''),',...
			'spm(''ChMod'',get(gco,''Value'')),',...
		'end'],...
					'Interruptible','off')

%-Utility buttons (second line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','',		'Position',[020 054 082 024].*WS,...
	'ToolTipString','<blank>',...
	'FontSize',FS2,			'CallBack','','Enable','off')

uicontrol(Fmenu,'Style','PopUp',...
	'String','Means...|Mean|adjMean|adjMean/fMRI',...
	'Position',[112 054 083 024].*WS,...
	'ToolTipString','image averaging utilities...',...
	'FontSize',FS2,			'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{	'spm_mean_ui',...
			'spm_adjmean_ui',...
			'spm_adjmean_fmri_ui'	})

uicontrol(Fmenu,'String','ImCalc',	'Position',[205 054 083 024].*WS,...
	'ToolTipString','image calculator',...
	'FontSize',FS2,			'CallBack','spm_imcalc_ui')

uicontrol(Fmenu,'String','HDR edit',	'Position',[298 054 082 024].*WS,...
	'ToolTipString','header editor',...
	'FontSize',FS2,			'CallBack','spm_header_edit')

%-Utility buttons (third line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Help',	'Position',[020 020 082 024].*WS,...
	'ToolTipString','launch SPMhelp browser',...
	'CallBack','spm_help',...
	'ForeGroundColor','g')

uicontrol(Fmenu,'Style','PopUp',...
	'String','Utils...|CD|Show SPM|Run mFile|SPMweb',...
	'ToolTipString','misc SPM utilities',...
	'Position',[112 020 083 024].*WS,...
	'FontSize',FS2,			'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{	[...
		'global CWD,',...
		'CWD=spm_get(-1,''*'',''Select new working directory'',pwd);',...
		'cd(CWD), clear CWD,',...
		'fprintf(''\nSPM working directory:\n\t%s\n\n>> '',pwd)'],...
			'spm(''Show'');',...
			'run(spm_get(1,''*.m'',''Select mFile to run''))',...
			'web(''http://www.fil.ion.ucl.ac.uk/spm'')' } )

uicontrol(Fmenu,'String','Defaults',	'Position',[205 020 083 024].*WS,...
	'ToolTipString','adjust default SPM behaviour for this session',...
	'FontSize',FS2,			'CallBack','spm_defaults_edit')

uicontrol(Fmenu,'String','Quit',	'Position',[298 020 082 024].*WS,...
	'ToolTipString','exit SPM',...
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
delete(spm_figure('FindWin','Interactive'))

%-Create SPM Interactive window
Rect = spm('WinSize','Interactive');
Finter = figure('IntegerHandle','off',...
	'Tag','Interactive',...
	'Name','','NumberTitle','off',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.7,...
	'MenuBar','none',...
	'DefaultUicontrolFontSize',2*round(10*min(spm('GetWinScale'))/2),...
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
% WS = spm('GetWinScale')
S   = get(0,'ScreenSize');
if all(S==1), error('Can''t open any graphics windows...'), end
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
% SPMver = spm('Ver',Mfile,ReDo,Cache)
if nargin<4, Cache=[]; else, Cache=varargin{4}; end
if isempty(Cache), Cache=1; end
if nargin<3, ReDo=[]; else, ReDo=varargin{3}; end
if isempty(ReDo), ReDo=0; end
if nargin<2, Mfile=''; else, Mfile=varargin{2}; end
if isempty(Mfile), Mfile='spm'; end
Fmenu = spm_figure('FindWin','Menu');

%-See if SPM Menu window exists - It's UserData should contain the version
%-----------------------------------------------------------------------
if ~ReDo & ~isempty(Fmenu)
	xSPM = get(Fmenu,'UserData');
	if isstruct(xSPM) & isfield(xSPM,'SPMver') & isfield(xSPM,'SPMc')
		varargout = {xSPM.SPMver,xSPM.SPMc};
		return
	end
end

%-Work version out from file
%-----------------------------------------------------------------------
SPMdir = spm('Dir',Mfile);
CFile  = [SPMdir,'/Contents.m'];
if exist(CFile)
	fid  = fopen(CFile,'r');
	str = setstr([fread(fid,80,'char')',setstr(10)]);
	fclose(fid);
	str(1:max(1,min(find(str~='%' & str~=' '))-1))=[];
	tmp = min(find(str==10|str==32));
	SPMver = str(1:tmp-1);
	if str(tmp)==32
		SPMc = str(tmp+1:tmp+min(find(str(tmp+1:end)==10)));
	else
		SPMc = '(c) The Wellcome Department of Cognitive Neurology';
	end
else
	SPMver = 'SPM';
	SPMc   = '(c) The Wellcome Department of Cognitive Neurology';
end

%-Store version info in UserData of SPM Menu window, if it exists
%-----------------------------------------------------------------------
if ~isempty(Fmenu) & Cache
	xSPM = get(Fmenu,'UserData');
	if ~isempty(xSPM) & ~isstruct(xSPM)
		warning('Fmenu UserData being overwritten!'), end
	xSPM.SPMver = SPMver;
	xSPM.SPMc   = SPMc;
	set(Fmenu,'UserData',xSPM);
end

varargout = {SPMver,SPMc};


case 'colour'
%=======================================================================
% spm('Colour')
%-----------------------------------------------------------------------
%-Developmental livery
varargout = {[0.7,1.0,0.7],'flourescent green'};
%-Alpha release livery
% varargout = {[0.9,0.9,0.5],'over-ripe banana'};
%-Beta release livery
% varargout = {[0.9 0.8 0.9],'blackcurrant purple'};
%-Distribution livery
% varargout = {[0.8 0.8 1.0],'vile violet'};


case 'getglobal'
%=======================================================================
% varargout = spm('GetGlobal',varargin)
wg = who('global');
for i=1:nargin-1
	if any(strcmp(wg,varargin{i+1}))
		eval(['global ',varargin{i+1},', tmp=',varargin{i+1},';'])
		varargout{i} = tmp;
	else
		varargout{i} = [];
	end
end

case 'isgcmdline'
%=======================================================================
% CmdLine = spm('isGCmdLine')
% if any(strcmp(who('global'),'CMDLINE')), global CMDLINE; else, CMDLINE=[]; end
CMDLINE = spm('GetGlobal','CMDLINE');
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


case 'popupcb'
%=======================================================================
% spm('PopUpCB',h)
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


case 'time'
%=======================================================================
% [timestr, date_vec] = spm('Time')
tmp = clock;
varargout = {sprintf('%02d:%02d:%02d - %02d/%02d/%4d',...
			tmp(4),tmp(5),floor(tmp(6)),tmp(3),tmp(2),tmp(1)),...
		tmp};


case 'pointer'
%=======================================================================
% spm('Pointer',Pointer)
if nargin<2, Pointer='Arrow'; else, Pointer=varargin{2}; end
set(get(0,'Children'),'Pointer',Pointer)


case 'fnbanner'
%=======================================================================
% SPMid = spm('FnBanner',Fn,FnV)
time = spm('time');
str  = spm('ver');
if nargin>=2, str = [str,': 'varargin{2}]; end
if nargin>=3, str = [str,' (v',varargin{3},')']; end
fprintf('\n%s',str)
fprintf('%c',' '*ones(1,72-length([str,time])))
fprintf('%s\n',time)
fprintf('%c','='*ones(1,72)),fprintf('\n')
varargout = {str};


case 'fnuisetup'
%=======================================================================
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup',Iname,bGX,CmdLine)
if nargin<4, CmdLine=spm('isGCmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, bGX=1; else, bGX=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end
if CmdLine
	Finter = spm_figure('FindWin','Interactive');
	if ~isempty(Finter), spm_figure('Clear',Finter), end
else
	Finter = spm_figure('GetWin','Interactive');
	spm_figure('Clear',Finter)
	str=spm('GetUser'); if ~isempty(Iname), str=[str,' - ',Iname]; end
	set(Finter,'Name',str)
end
if ~isempty(Iname), fprintf('%s:\n',Iname), end
if bGX
	Fgraph = spm_figure('GetWin','Graphics');
	spm_figure('Clear',Fgraph)
else
	Fgraph = spm_figure('FindWin','Graphics');
end
varargout = {Finter,Fgraph,CmdLine};	


case 'figname'
%=======================================================================
% F = spm('FigName',Iname,F,CmdLine)
if nargin<4, CmdLine=spm('isGCmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end

if ~isempty(Iname), fprintf('\t%s\n',Iname), end
if CmdLine, varargout={[]}; return, end
F = spm_figure('FindWin',F);
if ~isempty(F) & ~isempty(Iname)
	set(F,'Name',[spm('GetUser'),' - ',Iname])
end
varargout={F};


case 'show'
%=======================================================================
% Fs = spm('Show')
cF = get(0,'CurrentFigure');
Fs = get(0,'Children');
Fs = findobj(Fs,'flat','Visible','on');
for F=Fs', figure(F), end
set(0,'CurrentFigure',cF)
spm('FnBanner','GUI show');
varargout={Fs};


case 'clear'
%=======================================================================
% spm('Clear',Finter, Fgraph)
if nargin<3, Fgraph='Graphics'; else, Fgraph=varargin{3}; end
if nargin<2, Finter='Interactive'; else, Finter=varargin{2}; end
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)
spm('Pointer','Arrow')
spm_get('Initialise','reset');
clc, spm('FnBanner','GUI cleared');
fprintf('\n>> ');
%evalin('Base','clear')


case 'help'
%=======================================================================
% spm('Help',varargin)
if nargin>1, spm_help(varargin{2:end}), else, spm_help, end


otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end
