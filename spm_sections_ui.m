
% Display and analysis of regional effects
% FORMAT spm_sections_ui
%___________________________________________________________________________
%
% spm_sections_ui prompts for details of a SPM{Z} that is then
% displayed in the results window.  By moving (dragging) the pointer on
% the maximum intensity projection one can select a particular
% voxel (or planes passing though that voxel) for further display and
% analysis:
%
% 1) display activities at the point selected
% 2) plot the activities at the point selected
% 3) render the SPM(Z} on transverse slices at the level selected
% 4) render the SPM(Z} through 3 orthogonal sections selected
% 5) list and characterize all the maxima in a selected region
%
% The SPM{Z} is thresholded at a specified threshold and displayed
% on top of a background image.  This background image can be arbitrary
% but should have ORIGIN correctly specified in its header.  ORIGIN
% is the [i j k] voxel corresponding to [0 0 0] mm in XYZ [typically
% the anterior commissure in the space defined by the atlas of Talairach
% and Tournoux (1988)].
%
% The specification of the contrasts to use and uncorrected and corrected
% thresholds is the same as that described in spm_projections[_ui].  See
% SPM{Z} in the help facilitiy
%
%__________________________________________________________________________
% %W% %E%

% spm_sections_ui is not a function because many of the arguments and
% object handles required by CallBack strings and the spm_*.m routines
% that are invoked have to be in working memory

% Get directory, data, contrast and thresholds
%---------------------------------------------------------------------------
tmp  = spm_get(1,'.mat','select SPMt.mat for analysis','SPMt');
CWD  = strrep(tmp,'/SPMt.mat','');
K    = [];
FLIP = 1;
load([CWD,'/SPM' ])
load([CWD,'/XYZ' ])
load([CWD,'/SPMt'])
load([CWD,'/XA'  ])


con  = 0;
while (con < 1 | con > size(CONTRAST,1))
   con = spm_input(['which contrast[s] ? 1 - ' num2str(size(CONTRAST,1))],1);
end
U    = spm_input('threshold e.g. 2.8 or 0.001',2);
pV   = spm_input('corrected p value e.g. 0.05',3);


%---------------------------------------------------------------------------
set(2,'Pointer','watch'); figure(3); spm_clf

if U < 1
   U = spm_invNcdf(1 - U); end

% accommodate masking if required
%---------------------------------------------------------------------------
if length(con) > 1
	Q = all(SPMt(con,:) > U);
else
	Q = SPMt(con,:) > U;
end

% return if there are no voxels
%---------------------------------------------------------------------------
if sum(Q) == 0
	axis off
	text(0,0.8,CWD,'Fontsize',16,'FontWeight','Bold');
	text(0,0.7,'No voxels above this threshold');
	set(2,'Pointer','arrow'); return
end

t    = SPMt(con(1),Q);
XYZ  = XYZ(:,Q);
XA   = XA(:,Q);


% filter on P(nmax > k) and P(Zmax > u)
%---------------------------------------------------------------------------
A         = spm_clusters(XYZ,V([4 5 6]));
Q         = [];
for i     = 1:max(A)
	j = find(A == i);
	if ( spm_Pn(length(j),W,U,S) <= pV  ) | ( spm_Pz(W,max(t(j)),S) <= pV )
		Q = [Q j];
	end
end

% return if there are no voxels
%---------------------------------------------------------------------------
if sum(Q) == 0
	axis off
	text(0,0.8,CWD,'Fontsize',16,'FontWeight','Bold');
	text(0,0.7,'No voxels significant at this level');
	set(2,'Pointer','arrow'); return
end

t    = t(Q);
XYZ  = XYZ(:,Q);
XA   = XA(:,Q);


% display SPM{Z}
%---------------------------------------------------------------------------
figure(3); spm_clf;
set(3,'Units','normalized'); 
axes('Position',[0.24 0.54 0.62 0.42])
spm_mip(t,XYZ,V(1:6)); axis normal


% if 2-dimensional data
%---------------------------------------------------------------------------
if V(3) == 1
	set(gca,'Position',[0.24 0.54 0.52 0.36])
	X1   = text(1,1,'<','Color','r','Fontsize',16);
	X2   = X1;
	X3   = X1;

	% create 'clever' pointers that reset things when moved
	%-------------------------------------------------------------------
C1   = ['set(X1,''Units'',''Pixels'');',...
	'set(X1,''Position'',get(3,''CurrentPoint'') - AXES);',...
	'set(X1,''Units'',''Data'');',...
	'd   = get(X1,''Position'');',...
	'set(h1,''String'',sprintf(''%0.0f'',d(1)));'...
	'set(h2,''String'',sprintf(''%0.0f'',d(2)));'...
	'set(3,''WindowButtonUpFcn'','' '');',...
];
c    = ['d1   = eval(get(h1,''String''));',...
	'd2   = eval(get(h2,''String''));',...
	'set(X1,''Position'',[d1 d2 1]);'];

	C2    = C1;
	C3    = C1;
else

	% create line object = dot or selected voxel
	%------------------------------------------------------------------
	X1   = text(124,248,'<','Color','r','Fontsize',16);
	X2   = text(124,112,'<','Color','r','Fontsize',16);
	X3   = text(276,112,'<','Color','r','Fontsize',16);


	% create 'clever' pointers that reset things when moved
	%------------------------------------------------------------------
C1   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X1,''Units'',''Pixels'');',...
	'set(X1,''Position'',get(3,''CurrentPoint'') - AXES);',...
	'set(X1,''Units'',''Data'');',...
	'd   = get(X1,''Position'');',...
	'set(X2,''Position'',[d(1)     d2(2) 1]);',...
	'set(X3,''Position'',[(d3(1) + d(2) - d1(2)) d3(2) 1]);',...
	'set(h1,''String'',sprintf(''%0.0f'',(d(2) - 248)));'...
	'set(h2,''String'',sprintf(''%0.0f'',(d(1) - 124)));'...
	'set(3,''WindowButtonUpFcn'','' '');',...
];
C2   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X2,''Units'',''Pixels'');',...
	'set(X2,''Position'',get(3,''CurrentPoint'') - AXES);',...
	'set(X2,''Units'',''Data'');',...
	'd   = get(X2,''Position'');',...
	'set(X1,''Position'',[d(1)  d1(2) 1]);',...
	'set(X3,''Position'',[d3(1)  d(2) 1]);',...
	'set(h2,''String'',sprintf(''%0.0f'',(d(1) - 124)));'...
	'set(h3,''String'',sprintf(''%0.0f'',(112 - d(2))));'...
	'set(3,''WindowButtonUpFcn'','' '');',...
];
C3   = ['d1  = get(X1,''Position'');',...
	'd2  = get(X2,''Position'');',...
	'd3  = get(X3,''Position'');',...
	'set(X3,''Units'',''Pixels'');',...
	'set(X3,''Position'',get(3,''CurrentPoint'') - AXES);',...
	'set(X3,''Units'',''Data'');',...
	'd   = get(X3,''Position'');',...
	'set(X1,''Position'',[d1(1)  (d1(2) + d(1) - d3(1)) 1]);',...
	'set(X2,''Position'',[d2(1) d(2) 1]);',...
	'set(h1,''String'',sprintf(''%0.0f'',(d(1) - 276)));'...
	'set(h3,''String'',sprintf(''%0.0f'',(112 - d(2))));'...
	'set(3,''WindowButtonUpFcn'','' '');',...
];

c    = ['d1   = eval(get(h1,''String''));',...
	'd2   = eval(get(h2,''String''));',...
	'd3   = eval(get(h3,''String''));',...
	'set(X1,''Position'',[(124 + d2)  (248 + d1) 1]);',...
	'set(X2,''Position'',[(124 + d2)  (112 - d3) 1]);',...
	'set(X3,''Position'',[(276 + d1)  (112 - d3) 1]);'];
end

% reset axis attributes and get position in pixels
%---------------------------------------------------------------------------
set(gca,'Units','Pixels')
set(gcf,'Units','Pixels')
AXES = get(gca,'Position');
AXES = AXES(1:2);
set(gca,'Units','normalized')

% frames and text
%---------------------------------------------------------------------------
uicontrol(3,'Style','Frame','Position',[20 500 100 256]);
uicontrol(3,'Style','Text', 'Position',[30 720 16 16],'String','x');
uicontrol(3,'Style','Text', 'Position',[30 700 16 16],'String','y');
uicontrol(3,'Style','Text', 'Position',[30 680 16 16],'String','z');


% objects with callbacks that adjust the position of the dot
%---------------------------------------------------------------------------
set(X1,'ButtonDownFcn',['set(3,''WindowButtonUpFcn'',C1)']);
set(X2,'ButtonDownFcn',['set(3,''WindowButtonUpFcn'',C2)']);
set(X3,'ButtonDownFcn',['set(3,''WindowButtonUpFcn'',C3)']);

h1   = uicontrol(3,'Style','edit','String','1','Callback',c,...
         'Position',[60 720 40 16]);
h2   = uicontrol(3,'Style','edit','String','1','Callback',c,...
         'Position',[60 700 40 16]);
h3   = uicontrol(3,'Style','edit','String','1','Callback',c,...
         'Position',[60 680 40 16]);


% objects with calls to spm_*.m routines
%---------------------------------------------------------------------------
c    = ['L = [eval(get(h1,''String'')) eval(get(h2,''String'')) ',... 
	 'eval(get(h3,''String''))];'];

uicontrol(3,'Style','Pushbutton','String','table',   'Callback',...
	[c 'spm_table'   ],  'Position',[40 640 60 20],'Interruptible','yes');
uicontrol(3,'Style','Pushbutton','String','plot',    'Callback',...
	[c 'spm_plot'],      'Position',[40 610 60 20],'Interruptible','yes');
uicontrol(3,'Style','Pushbutton','String','slices','Callback',...
	[c 'spm_transverse'],'Position',[40 580 60 20],'Interruptible','yes');
uicontrol(3,'Style','Pushbutton','String','maxima','Callback',...
	[c 'spm_maxima'],    'Position',[40 520 60 20],'Interruptible','yes');

if V(3) > 1
	uicontrol(3,'Style','Pushbutton','String','sections',  'Callback',...
	[c 'spm_sections'],  'Position',[40 550 60 20],'Interruptible','yes');
end


% finished
%---------------------------------------------------------------------------
set(2,'Pointer','arrow')

















