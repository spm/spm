
% Display and analysis of SPM{F}
% FORMAT spm_spmF_ui
%___________________________________________________________________________
%
% spm_spmF_ui prompts for details of a SPM{F} that is then
% displayed in the results window.  By moving (dragging) the pointer on
% the maximum intensity projection one can select a particular
% voxel at which to plot the regional activities.
%
% The SPM{F} is thresholded at a specified (uncorrected) threshold
%
%__________________________________________________________________________
% %W% %E%

% spm_sections_ui is not a function because many of the arguments and
% object handles required by CallBack strings and the spm_*.m routines
% that are invoked have to be in working memory

% get data and contrast
%---------------------------------------------------------------------------
tmp = spm_get(1,'.mat','select SPMF.mat for analysis','SPMF');
CWD = strrep(tmp,'/SPMF.mat',''); % Get directory name
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMF'])

if isempty([B G]) | isempty([H C])
	df = [rank([H C B G]), df];
	if ~isempty(H), df=df-[1,0]; end %-See spm_spm!
else
	df   = [rank([H C]) df];
end
U    = spm_input('threshold e.g. 2.8 or 0.001',1);
if U < 1
	U = spm_invFcdf(1 - U,[df]); end

Q  = SPMF > U;

set(2,'Pointer','watch')

% display SPM{Z}
%---------------------------------------------------------------------------
figure(3); spm_clf;
set(3,'Units','normalized'); 
axes('Position',[0.24 0.54 0.62 0.42])
spm_mip(SPMF(Q),XYZ(:,Q),V(1:6)); axis normal
xlabel(sprintf('SPM{F}: p < %0.3f {uncorrected}',1 - spm_Fcdf(U,df)))
set(get(gca,'Xlabel'),'Visible','on')

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
[c 'spm_spmF_table'   ],    'Position',[40 640 60 20],'Interruptible','yes');
uicontrol(3,'Style','Pushbutton','String','plot',    'Callback',...
[c 'spm_spmF_plot'],        'Position',[40 610 60 20],'Interruptible','yes');

uicontrol(3,'Style','Pushbutton','String','hold',    'Callback',...
['subplot(2,1,2);hold on;'],'Position',[40 550 60 20],'Interruptible','yes');



% load adjusted data as aprelude to further analysis
%---------------------------------------------------------------------------
load([CWD,'/XA'])
load([CWD,'/BETA'])
set(2,'Pointer','arrow')

















