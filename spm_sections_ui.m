
% Display and analysis of regional effects
% FORMAT spm_sections_ui
%_______________________________________________________________________
%
% spm_sections_ui prompts for details of a SPM{Z} of SPM{F} that is then
% displayed in the results window.  By moving (dragging) the pointer on
% the maximum intensity projection one can select a particular
% voxel (or planes passing though that voxel) for further display and
% analysis:
%
% 1) plot the activities at the point selected
% 2) render the SPM(Z/F} on transverse slices at the level selected
% 3) render the SPM(Z/F} through 3 orthogonal sections selected
% 4) list and characterize all the maxima in a selected region
%
% The SPM{Z/F} is thresholded at a specified threshold and displayed
% on top of a background image.  This background image can be arbitrary
% but should have ORIGIN correctly specified in its header.  ORIGIN
% is the [i j k] voxel corresponding to [0 0 0] mm in XYZ [typically
% the anterior commissure in the space defined by the atlas of Talairach
% and Tournoux (1988)].
%
% The specification of the contrasts to use and the height and size
% thresholds are the same as that described in spm_projections_ui.%
%_______________________________________________________________________
% %W% %E%

% spm_sections_ui is not a function because many of the arguments and
% object handles required by CallBack strings and the spm_*.m routines
% that are invoked have to be in working memory

WS     = spm('GetWinScale');
Fgraph = spm_figure('FindWin','Graphics');
Finter = spm_figure('FindWin','Interactive');
spm_clf(Fgraph)
spm_clf(Finter)
global CWD

% which SPM
%-----------------------------------------------------------------------
SPMZ     = spm_input('which SPM',1,'b','SPM{Z}|SPM{F}',[1 0]);
SPMF     = ~SPMZ;

% get thresholded data, thresholds and parameters
%-----------------------------------------------------------------------
if SPMZ

	[t,XYZ,XA,U,k,s,w] = spm_projections_ui('Results');

elseif SPMF

	[t,XYZ,XA,U,k,s,w] = spm_projectionsF_ui('Results');

end
load([CWD,'/SPM'])
S        = s;
W        = w;

if SPMF
      df = Fdf; end

% display SPM{Z/F}
%=======================================================================
figure(Fgraph)
set(Fgraph,'Units','normalized'); 
axes('Position',[0.24 0.54 0.62 0.42])
spm_mip(t,XYZ,V(1:6)); axis normal


% if 2-dimensional data
%-----------------------------------------------------------------------
if V(3) == 1
	set(gca,'Position',[0.24 0.54 0.52 0.36])
	X1   = text(1,1,'<','Color','r','Fontsize',16);
	X1g  = text(1,1,'<','Color','g','Fontsize',12,'Visible','off');
	X2   = X1;
	X3   = X1;

	% create 'clever' pointers that reset things when moved
	%---------------------------------------------------------------
C1   = ['set(X1,''Units'',''Pixels'');',...
	'set(X1,''Position'',get(' num2str(Fgraph) ',''CurrentPoint'') - AXES);',...
	'set(X1,''Units'',''Data'');',...
	'd   = get(X1,''Position'');',...
	'set(h1,''String'',sprintf(''%0.0f'',d(1)));'...
	'set(h2,''String'',sprintf(''%0.0f'',d(2)));'...
	'set(hXstr,''String'',sprintf(''x = %0.0f'',d(1)));'...
	'set(hYstr,''String'',sprintf(''y = %0.0f'',d(2)));'...
	'set(' num2str(Fgraph) ',''WindowButtonUpFcn'','' '');',...
];
c    = ['d1   = eval(get(h1,''String''));',...
	'd2   = eval(get(h2,''String''));',...
	'set(X1,''Position'',[d1 d2 1]);'];

	C2    = C1;
	C3    = C1;
else
	% see spm_project.c for where these numbers come from
	%--------------------------------------------------------------
	P1 = 127-4;
	P2 = 182+182-91;
	P3 = 182-73;
	P4 = 218+91-4;

	% create line object = dot or selected voxel
	%--------------------------------------------------------------
	X1   = text(P1,P2,'<','Color','r','Fontsize',16);
	X2   = text(P1,P3,'<','Color','r','Fontsize',16);
	X3   = text(P4,P3,'<','Color','r','Fontsize',16);

	X1g  = text(P1,P2,'<','Color','g','Fontsize',12,'Visible','off');
	X2g  = text(P1,P3,'<','Color','g','Fontsize',12,'Visible','off');
	X3g  = text(P4,P3,'<','Color','g','Fontsize',12,'Visible','off');


	% create 'clever' pointers that reset things when moved
	%--------------------------------------------------------------

	% Transverse
	%--------------------------------------------------------------
C1   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X1,''Units'',''Pixels'');',...
	'set(X1,''Position'',get(' num2str(Fgraph) ',''CurrentPoint'') - AXES);',...
	'set(X1,''Units'',''Data'');',...
	'd   = get(X1,''Position'');',...
	'set(X2,''Position'',[d(1)     d2(2) 1]);',...
	'set(X3,''Position'',[(d3(1) + d(2) - d1(2)) d3(2) 1]);',...
	'set(h1,''String'',sprintf(''%0.0f'',(d(2) - ' num2str(P2) ')));'...
	'set(h2,''String'',sprintf(''%0.0f'',(d(1) - ' num2str(P1) ')));'...
	'set(hXstr,''String'',sprintf(''x = %0.0f'',(d(2) - ' num2str(P2) ')));'...
	'set(hYstr,''String'',sprintf(''y = %0.0f'',(d(1) - ' num2str(P1) ')));'...
	'set(' num2str(Fgraph) ',''WindowButtonUpFcn'','' '');',...
];

	% sagittal
	%--------------------------------------------------------------
C2   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X2,''Units'',''Pixels'');',...
	'set(X2,''Position'',get(' num2str(Fgraph) ',''CurrentPoint'') - AXES);',...
	'set(X2,''Units'',''Data'');',...
	'd   = get(X2,''Position'');',...
	'set(X1,''Position'',[d(1)  d1(2) 1]);',...
	'set(X3,''Position'',[d3(1)  d(2) 1]);',...
	'set(h2,''String'',sprintf(''%0.0f'',(d(1) - ' num2str(P1) ')));'...
	'set(h3,''String'',sprintf(''%0.0f'',(' num2str(P3) ' - d(2))));'...
	'set(hYstr,''String'',sprintf(''y = %0.0f'',(d(1) - ' num2str(P1) ')));'...
	'set(hZstr,''String'',sprintf(''z = %0.0f'',(' num2str(P3) ' - d(2))));'...
	'set(' num2str(Fgraph) ',''WindowButtonUpFcn'','' '');',...
];

	% coronal
	%--------------------------------------------------------------
C3   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X3,''Units'',''Pixels'');',...
	'set(X3,''Position'',get(' num2str(Fgraph) ',''CurrentPoint'') - AXES);',...
	'set(X3,''Units'',''Data'');',...
	'd   = get(X3,''Position'');',...
	'set(X1,''Position'',[d1(1)  (d1(2) + d(1) - d3(1)) 1]);',...
	'set(X2,''Position'',[d2(1) d(2) 1]);',...
	'set(h1,''String'',sprintf(''%0.0f'',(d(1) - ' num2str(P4) ')));'...
	'set(h3,''String'',sprintf(''%0.0f'',(' num2str(P3) ' - d(2))));'...
	'set(hXstr,''String'',sprintf(''x = %0.0f'',(d(1) - ' num2str(P4) ')));'...
	'set(hZstr,''String'',sprintf(''z = %0.0f'',(' num2str(P3) ' - d(2))));'...
	'set(' num2str(Fgraph) ',''WindowButtonUpFcn'','' '');',...
];

c    = ['d1   = eval(get(h1,''String''));',...
	'd2   = eval(get(h2,''String''));',...
	'd3   = eval(get(h3,''String''));',...
	'set(X1,''Position'',[(' num2str(P1) ' + d2)  (' num2str(P2) ' + d1) 1]);',...
	'set(X2,''Position'',[(' num2str(P1) ' + d2)  (' num2str(P3) ' - d3) 1]);',...
	'set(X3,''Position'',[(' num2str(P4) ' + d1)  (' num2str(P3) ' - d3) 1]);'];
end

% reset axis attributes and get position in pixels
%---------------------------------------------------------------------------
set(gca,'Units','Pixels')
set(gcf,'Units','Pixels')
AXES = get(gca,'Position');
AXES = AXES(1:2);
set(gca,'Units','normalized')

% print parameters
%---------------------------------------------------------------------------
hTextAxes = axes('Position',[0.625 0.54 0.3 0.15],'Visible','off',...
	'Units','Normalized');
text(0.0,1.0,spm('DirTrunc',CWD,24),'FontSize',8,'FontWeight','Bold')
text(0,0.8,sprintf('Height threshold {u} = %0.2f, p = %0.3f',U,1 - spm_Ncdf(U)),'FontSize',8)
text(0,0.7,...
	sprintf('Extent threshold {k} = %i voxels',k),...
	'FontSize',8)
hXstr = text(0,0.5,sprintf('x = 1'),...
	'FontSize',8);
hYstr = text(0,0.4,sprintf('y = 1'),...
	'FontSize',8);
hZstr = text(0,0.3,sprintf('z = 1'),...
	'FontSize',8);

% frames and text
%---------------------------------------------------------------------------
uicontrol(Fgraph,'Style','Frame','Position',[20 500 100 256].*WS);
uicontrol(Fgraph,'Style','Text', 'Position',[30 720 16 16].*WS,'String','x');
uicontrol(Fgraph,'Style','Text', 'Position',[30 700 16 16].*WS,'String','y');
uicontrol(Fgraph,'Style','Text', 'Position',[30 680 16 16].*WS,'String','z');


% objects with callbacks that adjust the position of the dot
%---------------------------------------------------------------------------
set(X1,'ButtonDownFcn',['set(' num2str(Fgraph) ',''WindowButtonUpFcn'',C1)']);
set(X2,'ButtonDownFcn',['set(' num2str(Fgraph) ',''WindowButtonUpFcn'',C2)']);
set(X3,'ButtonDownFcn',['set(' num2str(Fgraph) ',''WindowButtonUpFcn'',C3)']);

h1   = uicontrol(Fgraph,'Style','edit','String','1','Callback',c,...
         'Position',[60 720 40 16].*WS);
h2   = uicontrol(Fgraph,'Style','edit','String','1','Callback',c,...
         'Position',[60 700 40 16].*WS);
h3   = uicontrol(Fgraph,'Style','edit','String','1','Callback',c,...
         'Position',[60 680 40 16].*WS);


% objects with calls to spm_*.m routines
%---------------------------------------------------------------------------
c    = ['L = [eval(get(h1,''String'')) eval(get(h2,''String'')) ',... 
	 'eval(get(h3,''String''))];'];


uicontrol(Fgraph,'Style','Pushbutton','String','plot', 'Callback',...
	[c 'spm_graph'],  'Position',[40 610 60 20].*WS,'Interruptible','yes');
uicontrol(Fgraph,'Style','Pushbutton','String','slices','Callback',...
	[c 'spm_transverse'],'Position',[40 580 60 20].*WS,...
	'Interruptible', 'yes');
uicontrol(Fgraph,'Style','Pushbutton','String','maxima','Callback',...
	[c 'spm_maxima'],'Position',[40 520 60 20].*WS,'Interruptible','yes');

if V(3) > 1
	uicontrol(Fgraph,'Style','Pushbutton','String','sections',... 
	'Callback',[c 'spm_sections'],'Position',...
	[40 550 60 20].*WS,'Interruptible','yes');
end


% finished
%---------------------------------------------------------------------------
set(Finter,'Pointer','Arrow')

















