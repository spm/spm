function [R1,R2]=spm(Action,P2,P3)
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
% availabe to the scientific community.
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
% (equivalently `spm pet` and `spm mri` lead directly to the respective
% modality interfaces.
%
% Once the modlity is chosen, (and it can be toggled mid-session) the
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
% Arguments to this routine lead to various setup facilities, mainly of
% use to SPM power users and programmers. See the programmers FORMAT &
% help in the main body of spm.m
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
% FORMAT spm('SetWinDefaults')
% Sets defaults for figures, such as white backgrounds, a4 paper, gray
% colormaps etc...
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
% FORMAT SPMver=spm('Ver')
% Returns the current version of SPM, identified by the top line of the
% spm.man file in the directory containing the currently used spm.m
%
% FORMAT [c,cName]=spm('Colour')
% Returns the rgb tripls and a description for the current en-vogue SPm
% colour, the background colour for the Menu and Help windows.
%
% CmdLine = spm('isGCmdLine')
% Returns true if global CMDLINE exists and is itself true.
%
% FORMAT spm('SetCmdWinLabel',WinStripe,IconLabel)
% Sets the names on the headers and icons of Sun command tools.
% WinStripe defaults to a summary line identifying the user, host and
% MatLab version; IconLabel to 'MatLab'.
%
% FORMAT TDname = spm('DirTrunc',Dname,len)
% Truncates directory names to maximum length len characters,
% prepending '...' and keeping components intact.
%
%_______________________________________________________________________

%-Parameters
%-----------------------------------------------------------------------
Modalities = str2mat('PET','FMRI');

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Welcome'; end



if strcmp(lower(Action),lower('Welcome'))
%=======================================================================

%-Clear the Command window & delete all figures
%-----------------------------------------------------------------------
clc
close all

%-Print welcome banner
%-----------------------------------------------------------------------
spm('AsciiWelcome')

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
S = get(0,'ScreenSize');
F = figure('Color',[1 1 1]*.8,...
	'Name','',...
	'NumberTitle','off',...
	'Position',[S(3)/2-300,S(4)/2-140,500,280],...
	'Resize','off',...
	'Tag','Welcome',...
	'Pointer','Watch',...
	'Visible','off');
spm('SetWinDefaults')
spm('SetCmdWinLabel')

%-Frames and text
%-----------------------------------------------------------------------
axes('Position',[0 0 80/500 280/280],'Visible','Off')
text(0.5,0.5,'SPM',...
	'FontName','Times','FontSize',96,...
	'Rotation',90,...
	'VerticalAlignment','middle','HorizontalAlignment','center',...
	'Color',[1 1 1]*.6);

uicontrol(F,'Style','Frame','Position',[110 120 380 140]);
uicontrol(F,'Style','Frame','Position',[110 020 380 090]);

c = ['STATISTICAL PARAMETRIC MAPPING  - {',spm('Ver'),'}'];
uicontrol(F,'Style','Text', 'Position',[112 220 376 30],...
'String',c,'ForegroundColor',[1 0 0])

c = 'The Wellcome Department of Cognitive Neurology,';
uicontrol(F,'Style','Text', 'Position',[112 200 376 016],'String',c)
c = 'The Institute of Neurology';
uicontrol(F,'Style','Text', 'Position',[112 180 376 016],'String',c)
c = 'University College London';
uicontrol(F,'Style','Text', 'Position',[112 160 376 016],'String',c)

%-Objects with Callbacks - PET, fMRI or About SPM
%-----------------------------------------------------------------------
uicontrol(F,'String','PET and SPECT',...
	'Position',[140 066 150 030],...
	'CallBack',...
		'set(gcf,''Pointer'',''watch''), clear all, spm(''PET'')',...
	'Interruptible','yes',...
	'ForegroundColor',[0 1 1]);
uicontrol(F,'String','fMRI time-series',...
	'Position',[310 066 150 030],...
	'CallBack',...
		'set(gcf,''Pointer'',''watch''), clear all, spm(''FMRI'')',...
	'Interruptible','yes',...
	'ForegroundColor',[0 1 1]);
uicontrol(F,'String','About SPM',...
	'Position',[140 030 150 030],...
	'CallBack',[...
		'set(gcf,''Visible'',''off''),',...
		'spm_help(''spm.man'')'],...
	'Interruptible','yes',...
	'ForegroundColor','g');
uicontrol(F,'String','Quit',...
	'Position',[310 030 150 030],...
	'CallBack','close all,clear all,clc,fprintf(''Bye...\n\n>> '')',...
	'ForegroundColor','r');

set(F,'Pointer','Arrow','Visible','on')

return


elseif strcmp(lower(Action),lower('AsciiWelcome'))
%=======================================================================
disp( ' ___  ____  __  __                                                  ')
disp( '/ __)(  _ \(  \/  )  Statistical Parametric Mapping                 ')
disp( '\__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology ')
disp(['(___/(__)  (_/\/\_)  Version: ',spm('Ver')])
disp(' ')
% disp('  John Ashburner, Karl Friston, Andrew Holmes, Jean-Baptiste Poline')
fprintf('\n')
return


elseif strcmp(lower(Action),lower('PET')) | ...
	strcmp(lower(Action),lower('FMRI'))
%=======================================================================
% spm(Modality)

clc
spm('AsciiWelcome'),			fprintf('Initialising')

Modality = upper(Action);

close all,				fprintf('.')

%-Draw SPM windows
%-----------------------------------------------------------------------
Fmenu = spm('CreateMenuWin','off');	fprintf('.')
Finter = spm('CreateIntWin','off');	fprintf('.')
spm('SetWinDefaults')
Fgraph = spm_figure('Create','Graphics','Graphics','off');
					fprintf('.')

Fmotd = [spm('Dir'),'/spm_motd.man'];
if exist(Fmotd), spm_help('!Disp',Fmotd,'',Fgraph,spm('Ver'));
	else, spm_figure('WaterMark',Fgraph,spm('Ver')), end
					fprintf('.')

%-Setup for current modality
%-----------------------------------------------------------------------
spm('ChMod',Modality),			fprintf('.')

%-Reveal windows
%-----------------------------------------------------------------------
set([Fmenu,Finter,Fgraph],'Visible','on')
					fprintf('\n')

return




elseif strcmp(lower(Action),lower('CreateMenuWin'))
%=======================================================================
% Fmenu = spm('CreateMenuWin',Vis)
if nargin<2, Vis='on'; else, Vis=P2; end

%-Close any existing 'Menu' 'Tag'ged windows
close(spm_figure('FindWin','Menu'))

%-Get size and scalings and create Menu window
%-----------------------------------------------------------------------
A = spm('GetWinScale');
Rect = spm('WinSize','Menu','raw').*A;

Fmenu = figure('Tag','Menu',...
	'Name',spm('Ver'),...
	'Color',[1 1 1]*.8,...
	'Position',Rect,...
	'NumberTitle','off',...
	'Resize','off',...
	'Visible','off');

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 145 380 295].*A,'Tag','Empty');

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 320 360 110].*A,'Tag','Empty');	
uicontrol(Fmenu,'Style','Text','String','Spatial',...
	'Position',[100 400 200 20].*A,'Tag','Empty',...
	'ForegroundColor','w')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 235 360 075].*A,'Tag','Empty');
	uicontrol(Fmenu,'Style','Text','String','Analysis',...
	'Position',[100 285 200 20].*A,'Tag','Empty',...
	'ForegroundColor','w')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 155 360 070].*A,'Tag','Empty');
uicontrol(Fmenu,'Style','Text','String','Results',...
	'Position',[100 200 200 20].*A,'Tag','Empty',...
	'ForegroundColor','w' )

uicontrol(Fmenu,'Style','Text',...
	'String','SPM for PET/SPECT',...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'HorizontalAlignment','center',...
	'Position',[020 125 360 020].*A,...
	'Tag','PET','Visible','off')
uicontrol(Fmenu,'Style','Text',...
	'String','SPM for functional MRI',...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'HorizontalAlignment','center',...
	'Position',[020 125 360 020].*A,...
	'Tag','FMRI','Visible','off')

uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 010 380 112].*A,'Tag','Empty');

%-Objects with Callbacks - main spm_*_ui.m routines
%-----------------------------------------------------------------------

%-Spatial
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Realign',	'Position',[040 370 080 30].*A,...
	'CallBack','spm_realign',	'Interruptible','yes');

uicontrol(Fmenu,'String','Normalize',	'Position',[150 350 100 30].*A,...
	'CallBack','spm_sn3d',	'Interruptible','yes');

uicontrol(Fmenu,'String','Smooth',	'Position',[280 370 080 30].*A,...
	'CallBack','spm_smooth_ui',	'Interruptible','yes');

uicontrol(Fmenu,'String','Coregister',	'Position',[040 330 080 30].*A,...
	'CallBack','spm_coregister',	'Interruptible','yes');

uicontrol(Fmenu,'String','Segment',	'Position',[280 330 080 30].*A,...
	'CallBack','spm_segment',	'Interruptible','yes');

%-Statistics
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Statistics',	'Position',[040 245 140 30].*A,...
	'CallBack','spm_spm_ui',	'Interruptible','yes',...
	'Visible','off',		'Tag','PET');

uicontrol(Fmenu,'String','Statistics',	'Position',[040 245 140 30].*A,...
	'CallBack','spm_fmri_spm_ui',   'Interruptible','yes',...
	'Visible','off',		'Tag','FMRI');

uicontrol(Fmenu,'String','Eigenimages',	'Position',[220 245 140 30].*A,...
	'CallBack','spm_svd_ui',	'Interruptible','yes');

%-Results
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','SPM{F}',	'Position',[045 165 070 30].*A,...
	'CallBack','spm_projectionsF_ui','Interruptible','yes');

uicontrol(Fmenu,'String','Results',	'Position',[165 165 070 30].*A,...
	'CallBack','spm_results_ui',	'Interruptible','yes');

uicontrol(Fmenu,'String','SPM{Z}',	'Position',[285 165 070 30].*A,...
	'CallBack','spm_projections_ui','Interruptible','yes');

%-Utility buttons (first line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Analyze',	'Position',[020 088 082 024].*A,...
	'CallBack','!analyze',		'Interruptible','yes');

uicontrol(Fmenu,'String','Display',	'Position',[112 088 083 024].*A,...
	'CallBack','spm_image',		'Interruptible','yes');

uicontrol(Fmenu,'String','Render',	'Position',[205 088 083 024].*A,...
	'CallBack','spm_render',	'Interruptible','yes');

uicontrol(Fmenu,'Style','PopUp','String',Modalities,...
	'Tag','Modality',		'Position',[298 088 082 024].*A,...
	'CallBack',...
	[	'if get(gco,''Value'')~=get(gco,''UserData''),',...
			'spm(''ChMod'',get(gco,''Value'')),',...
		'end'],...
					'Interruptible','yes');

%-Utility buttons (second line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','GhostView',	'Position',[020 054 082 024].*A,...
	'CallBack',[...
		'unix([''ghostview '',spm_get(1,''.ps'',''Select ',...
		'PostScript file to view''),'' &'']);',...
			],		'Interruptible','yes');

uicontrol(Fmenu,'String','CD',		'Position',[112 054 083 024].*A,...
	'CallBack',[...
		'cd(spm_get(-1,[],''Select new working directory'')),',...
		'fprintf(''\nSPM working directory:\n\t%s\n\n>> '',pwd)',...
			],		'Interruptible','yes');

uicontrol(Fmenu,'String','Mean',	'Position',[205 054 083 024].*A,...
	'CallBack','spm_average',	'Interruptible','yes');

uicontrol(Fmenu,'String','ImCalc',	'Position',[298 054 082 024].*A,...
	'CallBack','spm_image_funks',	'Interruptible','yes');

%-Utility buttons (third line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Help',	'Position',[020 020 082 024].*A,...
	'CallBack','spm_help',		'Interruptible','yes',...
	'ForeGroundColor','g');

uicontrol(Fmenu,'String','Defaults',	'Position',[112 020 083 024].*A,...
	'CallBack','spm_defaults_edit',	'Interruptible','yes');

User = getenv('USER');
uicontrol(Fmenu,'String',['<',User,'>'],'Position',[205 020 083 024].*A,...
					'Interruptible','yes',...
	'CallBack',...
		['if exist(''' User ''');',...
		 User,';else;spm_help(''UserButton''); end']);

uicontrol(Fmenu,'String','Quit',	'Position',[298 020 082 024].*A,...
	'ForeGroundColor','r',...
	'CallBack','close all,clear all,clc,fprintf(''Bye...\n\n>> '')');
%-----------------------------------------------------------------------
% set(Fmenu,'Visible',Vis)
R1 = Fmenu;
return



elseif strcmp(lower(Action),lower('CreateIntWin'))
%=======================================================================
% Finter = spm('CreateIntWin',Vis)
if nargin<2, Vis='on'; else, Vis=P2; end

%-Close any existing 'Interactive' 'Tag'ged windows
close(spm_figure('FindWin','Interactive'))

%-Create SPM Interactive window
Rect = spm('WinSize','Interactive');
Finter = figure('Tag','Interactive',...
	'Name','',...
	'Color',[1 1 1]*.7,...
	'Position',Rect,...
	'NumberTitle','off',...
	'Resize','off',...
	'Visible',Vis);
R1 = Finter;
return


elseif strcmp(lower(Action),lower('ChMod'))
%=======================================================================
% spm('ChMod',Modality)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality = ''; else, Modality = P2; end
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
	if (w < 256) warndlg(['To increase; quit SPM & MatLab and type ',...
		'''unlimit'' and restart'],...
		sprintf('WARNING: file descriptors = %d',w))
	end
end
return




elseif strcmp(lower(Action),lower('SetWinDefaults'))
%=======================================================================
% spm('SetWinDefaults')
whitebg(0,'w')
set(0,'DefaultFigureColormap',gray);
set(0,'DefaultFigurePaperType','a4letter')

return

elseif strcmp(lower(Action),lower('defaults'))
%=======================================================================
% spm('defaults',Modality)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=P2; end
Modality = spm('CheckModality',Modality);

%-Set MODALITY
%-----------------------------------------------------------------------
global MODALITY
MODALITY = Modality;

%-Set global defaults (global variables)
%-----------------------------------------------------------------------
global SWD TWD DESCRIP UFp
global CWD

SWD	 = spm('Dir');					% SPM directory
CWD	 = pwd;						% Working directory
TWD	 = getenv('SPMTMP');				% Temporary directory
if isempty(TWD)
	TWD = '/tmp';
end

% Load system defaults
%-----------------------------------------------------------------------
global PET_DIM PET_VOX PET_TYPE PET_SCALE PET_OFFSET PET_ORIGIN PET_DESCRIP
global fMRI_DIM fMRI_VOX fMRI_TYPE fMRI_SCALE fMRI_OFFSET fMRI_ORIGIN fMRI_DESCRIP

spm_defaults;

UFp	 = 0.05;					% Upper tail F-prob

%-Set Modality specific default (global variables)
%-----------------------------------------------------------------------
global DIM VOX SCALE TYPE OFFSET ORIGIN
if strcmp(Modality,'PET')
	DIM	= PET_DIM;				% Dimensions [x y z]
	VOX	= PET_VOX;				% Voxel size [x y z]
	SCALE	= PET_SCALE;				% Scaling coeficient
	TYPE	= PET_TYPE;				% Data type
	OFFSET	= PET_OFFSET;		    		% Offset in bytes
	ORIGIN	= PET_ORIGIN;				% Origin in voxels
elseif strcmp(Modality,'FMRI')
	DIM	= fMRI_DIM;				% Dimensions {x y z]
	VOX	= fMRI_VOX;				% Voxel size {x y z]
	SCALE	= fMRI_SCALE;				% Scaling coeficient
	TYPE	= fMRI_TYPE;				% Data type
	OFFSET	= fMRI_OFFSET;	   			% Offset in bytes
	ORIGIN	= fMRI_ORIGIN;				% Origin in voxels
elseif strcmp(Modality,'UNKNOWN')
else
	error('Illegal Modality')
end

return

elseif strcmp(lower(Action),lower('CheckModality'))
%=======================================================================
% [Modality,ModNum] = spm('CheckModality',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=upper(P2); end
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
R1 = upper(Modality);
R2 = ModNum;
return


elseif strcmp(lower(Action),lower('GetWinScale'))
%=======================================================================
% spm('GetWinScale')
S   = get(0,'ScreenSize');
R1  = [S(3)/1152 S(4)/900 S(3)/1152 S(4)/900];

return


elseif strcmp(lower(Action),lower('WinSize'))
%=======================================================================
% Rect = spm('WinSize',Win,raw)
if nargin<3, raw=0; else, raw=1; end
if nargin<2, Win=''; else, Win=P2; end

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
R1 = Rect;
return


elseif strcmp(lower(Action),lower('Dir'))
%=======================================================================
% spm('Dir',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='spm'; else, Mfile=P2; end
tmp    = which(Mfile);
if ~isstr(tmp)
	error(['Can''t find ',Mfile,' on MATLABPATH']); end
SPMdir = strrep(tmp,['/',Mfile,'.m'],'');
R1     = SPMdir;
return

elseif strcmp(lower(Action),lower('Ver'))
%=======================================================================
% spm('Ver',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='spm'; else, Mfile=P2; end
SPMdir = spm('Dir',Mfile);
CFile  = [SPMdir,'/spm.man'];
if exist(CFile)
	fid  = fopen(CFile,'r');
	SPMver = setstr([fread(fid,80,'char')',setstr(10)]);
	fclose(fid);
	SPMver(1:max(1,min(find(SPMver~='%' & SPMver~=' '))-1))=[];
	SPMver = SPMver(1:min(find(SPMver==10))-1);
else
	SPMver = 'SPM';
end
R1 = SPMver;
return

elseif strcmp(lower(Action),lower('Colour'))
%=======================================================================
% spm('Colour')
%-----------------------------------------------------------------------
%-Developmental livery
% R1 = [0.7,1.0,0.7];
% R2 = 'Lime Green';
%-Distribution livery
R1 = [0.8 0.8 1.0];
R2 = 'Diluted Blackcurrent Purple';
return


elseif strcmp(lower(Action),lower('isGCmdLine'))
%=======================================================================
% CmdLine = spm('isGCmdLine')
global CMDLINE
if isempty(CMDLINE), R1 = 0; else, R1 = CMDLINE; end
return


elseif strcmp(lower(Action),lower('SetCmdWinLabel'))
%=======================================================================
% spm('SetCmdWinLabel',WinStripe,IconLabel)
%-----------------------------------------------------------------------

%-Only label Sun command tools
%-----------------------------------------------------------------------
Term        = getenv('TERM');
if ~strcmp(Term,'sun-cmd'), return, end

%-Work out label text
%-----------------------------------------------------------------------
User        = getenv('USER');
[null,Host] = unix('echo `hostname` | sed -e ''s/\..*$//''');
Host        = Host(1:length(Host)-1);

if nargin<3, IconLabel = 'MatLab'; end
if nargin<2, WinStripe = [User,' - ',Host,' : MatLab ',version]; end

%-Set window stripe
%-----------------------------------------------------------------------
disp([']l' WinStripe '\]L' IconLabel '\'])
return

elseif strcmp(lower(Action),lower('DirTrunc'))
%=======================================================================
% TDname = spm('DirTrunc',Dname,len)
%-----------------------------------------------------------------------
if nargin < 3, len = 50; else, len = P3; end
if nargin < 2, Dname = spm('Dir'); else, Dname = P2; end
if length(Dname) > len
	lDname = length(Dname);
	tmp    = min(find(Dname(lDname-len:lDname)=='/')) + lDname-len -1;
	if isempty(tmp), tmp=max(find(Dname=='/')); end
	TDname = ['...',Dname(tmp:lDname)];
else
	TDname=Dname;
end

R1    = TDname;
return


else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
