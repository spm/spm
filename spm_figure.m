function spm_figure
% Creates the menu bar for graphics window [usually figure(3)]
% FORMAT spm_figure
%____________________________________________________________________________
%
% spm_figure creates user interface objects in the 'results window'
% that faciliate interactive editing of graphic output prior to
% printing (e.g. selection of color maps, deleting, moving and
% editing graphics objects or adding text)
% Objects with the attribute Tag = 'NoDelete' are exempt from deletion
% when spm_clf is used
%
% see also: spm_clf.m, spm_graphics_cut.m, spm_graphics_drag.m, 
% spm_graphics_size.m and spm_graphics_text.m
%
%__________________________________________________________________________
% %W% %E%


% get position parameters
%----------------------------------------------------------------------------
figure(3); clf; set(3,'Units','Pixels');
P  = get(3,'Position');
P  = P(3:4);					% figure dimensions {pixels}
x  = 50;					% uicontrol object width
y  = 20;					% uicontrol object height
z  = 0;						% initial x position

% create uicontrol objects
%----------------------------------------------------------------------------
uicontrol(3,'Style', 'Frame', 'Position',[0 (P(2) - 25) P(1) 30],...
	'Tag','NoDelete'); P = P(2) - y;

uicontrol(3,'String','Print' ,'Position',[z  P  60 y],'CallBack',...
	'spm_print',...
	'Tag','NoDelete','Foregroundcolor',[1 0 0]); z = z + 70;

uicontrol(3,'String','Clear' ,'Position',[z  P  60 y],'CallBack',...
	 ['figure(2); clf; figure(3); clf; spm_figure;',...
	  'set(2,''Name'','' '',''Pointer'',''Arrow'');',...
	  'set(3,''Pointer'',''Arrow'');'],...
        'Tag','NoDelete','Foregroundcolor',[1 0 0]); z = z + 70;

uicontrol(3,'String','gray'  ,'Position',[z  P  x  y],'CallBack',...
	'colormap(gray(64))',         'Tag','NoDelete'); z = z + x;
uicontrol(3,'String','hot'   ,'Position',[z  P  x  y],'CallBack',...
	'colormap(hot(64))' ,         'Tag','NoDelete'); z = z + x;
uicontrol(3,'String','split' ,'Position',[z  P  x  y],'CallBack',...
	'load Split; colormap(split)','Tag','NoDelete'); z = z + x;
uicontrol(3,'String','invert','Position',[z P  x y],'CallBack',...
        'colormap(flipud(colormap))', 'Tag','NoDelete'); z = z + x + 30;

uicontrol(3,'String','cut',   'Position',[z P  x  y],'CallBack',...
        'spm_graphics_cut',   'Tag','NoDelete'); z = z + x;
uicontrol(3,'String','move',  'Position',[z P  x  y],'CallBack',...
        'spm_graphics_drag',  'Tag','NoDelete'); z = z + x;
uicontrol(3,'String','size',  'Position',[z P  x  y],'CallBack',...
        'spm_graphics_size',  'Tag','NoDelete'); z = z + x;
uicontrol(3,'String','text',  'Position',[z P  x  y],'CallBack',...
        'spm_graphics_text',  'Tag','NoDelete'); z = z + x;

%----------------------------------------------------------------------------
colormap gray
spm_clf
