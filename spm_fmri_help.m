function spm_fmri_help
% Interactive help and manual facility
% FORMAT spm_fmri_help
%______________________________________________________________________________
%
% Provides a layout of the SPM components detailing their relation
% to each other.  Help is available on each of these components and
% any of the routines used to implement them.
% This window and its children are configured at the beginning of
% the session and is made accessible using the figure's 'Visible'
% attribute
%
%__________________________________________________________________________
% %W% %E%
 


% position figures and set PaperPosition for printing
%----------------------------------------------------------------------------
S     = get(0,'ScreenSize');
A     = diag([S(3)/1152 S(4)/900 S(3)/1152 S(4)/900]);
S1    = [108 429 400 445]*A;
F     = figure('Color',[0.8 0.8 1],'Name','SPM95 {fMRI} Help facility',...
	'NumberTitle','off','Position',S1,'Resize','off','Visible','off');

COLOR = [0 0 1];


% objects with Callbacks - spm_sub(.) routines
%----------------------------------------------------------------------------
c  = 'spm_sub(''spm_bib.man'');'  ;
uicontrol(F,'String','About SPM',  'Position',[010 395 90 30]*A,'CallBack',c,'ForegroundColor',COLOR);

c  = 'spm_sub(''spm_format.man'');'     ;
uicontrol(F,'String','Data format','Position',[110 395 90 30]*A,'CallBack',c,'ForegroundColor',COLOR);

c  = 'spm_sub(''spm_fmri.man'');'    ;
uicontrol(F,'String','Overview',   'Position',[210 395 90 30]*A,'CallBack',c,'ForegroundColor',[1 0 0]);

c  = 'spm_sub(''spm_ui.man'');'         ;
uicontrol(F,'String','Variables',  'Position',[310 395 86 30]*A,'CallBack',c,'ForegroundColor',COLOR);

c  = 'spm_sub(''spm_results.man'');'    ;
uicontrol(F,'String','Graphics',   'Position',[160 090 070 24]*A,'CallBack',c,'ForegroundColor',COLOR);


c  = 'spm_sub(''spm_realign.man'');'    ;
uicontrol(F,'String','Realign',    'Position',[040 350 080 30]*A,...
	'CallBack',c,    		'Interruptible','yes');

c  = 'spm_sub(''spm_sn.man'');'         ;
uicontrol(F,'String','Normalize',  'Position',[150 350 100 30]*A,...
	'CallBack',c,         		'Interruptible','yes');

c  = 'spm_sub(''spm_smooth.man'');'     ;
uicontrol(F,'String','Smooth',     'Position',[280 350 080 30]*A,...
	'CallBack',c,     		'Interruptible','yes');

c  = 'spm_sub(''spm_spm.man'');'        ;
uicontrol(F,'String','Statistics', 'Position',[130 270 140 30]*A,...
	'CallBack',c,        	   'Interruptible','yes');

c  = 'spm_sub(''spm_svd.man'');'        ;
uicontrol(F,'String','Eigenimages','Position',[130 220 140 30]*A,...
	'CallBack',c,        	   'Interruptible','yes');

c  = 'spm_sub(''spm_F.man'');'         ;
uicontrol(F,'String','SPM{F}',     'Position',[045 130 070 30]*A,...
	'CallBack',c,          	   'Interruptible','yes');

c  = 'spm_sub(''spm_projections.man'');';
uicontrol(F,'String','SPM{Z}',     'Position',[165 130 070 30]*A,...
	'CallBack',c,		   'Interruptible','yes');

c  = 'spm_sub(''spm_sections.man'');'   ;
uicontrol(F,'String','Results',    'Position',[285 130 070 30]*A,...
	'CallBack',c,   	   'Interruptible','yes');

uicontrol(F,'String','Analyze',    'Position',[020 56 80 24]*A,...
	'CallBack',' ',         	   'Interruptible','yes');

c  = 'spm_sub(''spm_image.man'');'      ;
uicontrol(F,'String','display',    'Position',[110 56 80 24]*A,...
	'CallBack',c,         	   'Interruptible','yes');

c  = 'spm_sub(''spm_average.m'');'      ;
uicontrol(F,'String','mean',       'Position',[200 56 80 24]*A,...
	'CallBack',c,              'Interruptible','yes');

c  = 'spm_sub(''spm_help.man'');'       ;
uicontrol(F,'String','help',       'Position',[020 20 80 24]*A,...
	'CallBack',c,      	   'Interruptible','yes');

c  = 'spm_sub(''spm_defaults.man'');'   ;
uicontrol(F,'String','defaults',   'Position',[110 20 80 24]*A,...
	'CallBack',c,      	   'Interruptible','yes');

c  = 'spm_sub(''spm_button.man'');'   ;
uicontrol(F,'String',getenv('USER'),'Position',[200 20 80 24]*A,...
	'CallBack',c,              'Interruptible','yes');

c  = 'set(gcf,''Visible'',''off''); figure(2); clf; figure(3); spm_clf';
uicontrol(F,'String','Quit',       'Position',[300 20 80 24]*A,...
	'CallBack',c);


set(F,'Tag','HELP')
