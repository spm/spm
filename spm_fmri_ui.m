function spm_fmri_ui
% SPM user interface
% FORMAT spm_fmri_ui
%____________________________________________________________________________
% spm_fmri_ui sets the default value of global variables and configures three
% Figures (children of root = 0):
%
% Figure 1 - pushbutton menu calling SPM routines
% Figure 2 - a panel for interactive specification of variables
% Figure 3 - the results window with various edit and print facilities
%
% These figures or windows do not change during an SPM session and provide
% a constant visual environment in which data analysis is implemented. The
% layout has been designed to be simple and at the same time show all the
% facilities that are available.
% Each user control object (pushbutton) has an associated 'Callback' or
% string that is evaluated on selection. These Callbacks are SPM routines
%
%__________________________________________________________________________
% %W% Karl Friston %E%



% Delete exisiting variables and figures
%----------------------------------------------------------------------------
clear; clc
close  all
clear global


% Set default values for gobal variables
%----------------------------------------------------------------------------
global SWD CWD PRINTSTR LOGFILE CMDLINE GRID DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP TWD UFp

SWD      = spm_get_path;				% SPM directory
CWD      = pwd;						% working directory
TWD      = getenv('SPMTMP');				% temporary directory
           if (TWD == []); TWD = '/tmp';end
PRINTSTR = ['print -dpsc -append ' CWD '/spm.ps'];    	% print string
LOGFILE  = '';						% SPM log file
if CMDLINE == [];                                       % command line input
        CMDLINE  = 0; end                               %   don't reset value
GRID     = 0.6;			    			% grid value

DIM	 = [128 64 10];					% dimensions {x y z]
VOX      = [3 3 5];					% voxel size {x y z]
SCALE    = 1;						% scaling coeficient
TYPE     = 4;						% data type
OFFSET   = 0;			    			% offset in bytes
ORIGIN   = [0 0 0];					% origin in voxels
DESCRIP	 = 'spm compatible';				% description

UFp      = 0.05;					% Upper tail F-prob

USER     = getenv('USER');

% position figures and set PaperPosition for printing
%----------------------------------------------------------------------------
S    = get(0,'ScreenSize');
A    = diag([S(3)/1152 S(4)/900 S(3)/1152 S(4)/900]);
S1   = [108 429 400 445]*A;
S2   = [108 008 400 395]*A;
S3   = [515 008 600 865]*A;

whitebg(0,'w')

figure('Color',[1 1 1]*.8,'Name','SPM for functional MRI {SPM95}',...
	'NumberTitle','off','Position',S1,'Resize','off','Visible','off');
figure('Color',[1 1 1]*.7,'Name','',...
	'NumberTitle','off','Position',S2,'Resize','off','Visible','off');
figure('Name','Results',...
	'NumberTitle','off','Position',S3,'Resize','off','Visible','off',...
	'PaperPosition',[.75 1.5 7 9.5])


% set up figure 4 - a help window that is usually invisible
%----------------------------------------------------------------------------
spm_fmri_help


% main menu (in Figure 1) - frames and text
%----------------------------------------------------------------------------
uicontrol(1,'Style','Frame','Position',[10 010 380 080]*A,...
'BackgroundColor',[.8 .8 1]);
uicontrol(1,'Style','Frame','Position',[10 100 380 340]*A,...
'BackgroundColor',[.8 .8 1]);
uicontrol(1,'Style','Frame','Position',[20 110 360 080]*A);
uicontrol(1,'Style','Frame','Position',[20 200 360 140]*A);
uicontrol(1,'Style','Frame','Position',[20 350 360 080]*A);

uicontrol(1,'Style','Text', 'Position',[100 405 200 20]*A,'String','spatial','ForegroundColor',[1 0 0] )
uicontrol(1,'Style','Text', 'Position',[100 315 200 20]*A,'String','analysis','ForegroundColor',[1 0 0])
uicontrol(1,'Style','Text', 'Position',[100 165 200 20]*A,'String','results','ForegroundColor',[1 0 0] )


% objects with Callbacks - main spm_*_ui.m routines
%----------------------------------------------------------------------------
uicontrol(1,'String','Realign',         'Position',[040 370 080 30]*A,...
	'CallBack','spm_fmri_realign_ui','Interruptible','yes');

uicontrol(1,'String','Normalize',       'Position',[150 370 100 30]*A,...
	'CallBack','spm_mri_sn_ui',     'Interruptible','yes');

uicontrol(1,'String','Smooth',          'Position',[280 370 080 30]*A,...
	'CallBack','spm_smooth_ui',     'Interruptible','yes');

uicontrol(1,'String','Statistics',      'Position',[130 270 140 30]*A,...
	'CallBack','spm_fmri_spm_ui',   'Interruptible','yes');

uicontrol(1,'String','Eigenimages',     'Position',[130 220 140 30]*A,...
	'CallBack','spm_svd_ui',        'Interruptible','yes');

uicontrol(1,'String','SPM{F}',          'Position',[045 130 070 30]*A,...
	'CallBack','spm_spmF',          'Interruptible','yes');

uicontrol(1,'String','SPM{Z}',          'Position',[165 130 070 30]*A,...
	'CallBack','spm_projections_ui','Interruptible','yes');

uicontrol(1,'String','Results',         'Position',[285 130 070 30]*A,...
	'CallBack','spm_sections_ui',   'Interruptible','yes');

uicontrol(1,'String','Analyze',         'Position',[020 56 80 24]*A,...
	'CallBack','!analyze',          'Interruptible','yes');

uicontrol(1,'String','display',         'Position',[110 56 80 24]*A,...
	'CallBack','spm_image',         'Interruptible','yes');

uicontrol(1,'String','mean',            'Position',[200 56 80 24]*A,...
	'CallBack','spm_average',       'Interruptible','yes');

uicontrol(1,'String','defaults',        'Position',[110 20 80 24]*A,...
	'CallBack','spm_defaults',      'Interruptible','yes');

uicontrol(1,'String','Quit',            'Position',[300 20 80 24]*A,...
	'CallBack','close all, clear all, clc');

c  = ['set(findobj(get(0,''Children''),''flat'',''Tag'',''HELP''),''Visible'', ''on'')'];
uicontrol(1,'String','help',            'Position',[020 20 80 24]*A,...
	'CallBack',c);

c  = ['if exist(''' USER ''');' USER '; else;spm_help_disp(''spm_BUTTON.m''); end'];
uicontrol(1,'String',USER,              'Position',[200 20 80 24]*A,...
	'CallBack',c,                   'Interruptible','yes');


% Tag figures and tag as 'Empty', objects in main menu (figure 1)
%----------------------------------------------------------------------------
set(1,'Tag','Menu')
set(2,'Tag','Interactive')
set(3,'Tag','Graphics')
set(get(1,'Children'),'Tag','Empty')

% reveal windows and configure graphics bar
%----------------------------------------------------------------------------
set(1,'Visible','on')
set(2,'Visible','on')
set(3,'Visible','on')
spm_figure


% check maximium number of files open (descriptors)
%----------------------------------------------------------------------------
[s,w] = unix('limit');
d     = findstr(w,'descriptors');
w     = eval(w((d + 11):length(w)));
if w < 256
  	c = 'To increase; quit MATLAB and type ''unlimit''';
	d = sprintf('WARNING descriptors = %0.0i',w);
	dialog('Style','warning','TextString',...
        c,'Name',d,'Position',[409 414 288 76]);
end

