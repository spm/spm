function varargout = spm_DesRep(varargin)
% Design reporting utilities
% FORMAT varargout = spm_DesRep(action,varargin)
%_______________________________________________________________________
%
% spm_DesRep is a multi-function function providing a suite of utility
% functions for various graphical reports on a given experimental
% design, embodied in a design matrix and other associated data
% structures.
%
%                           ----------------
%
% FORMAT spm_DesRep('Files&Factors',P,I,xC,sF,xs)
% Produces multi-page listing of files, factor indices, and covariates.
% P   - nx1 CellStr of filenames (i.e. {V.fname}')
% I   - nx4 matrix of factor indices
% xC  - Covariate structure array (see spm_spm_ui.m for definitions)
%       ('rc' & 'cname' fields used)
% sF  - 1x4 CellStr of factor names (i.e. D.sF of spm_spm_ui)
% xs  - [optional] structure of extra strings containing descriptive
%       information which is printed out after the files & variables listing.
%       The field names are used as sub-headings, the field values (which
%       must be strings or CellStr) printed alongside.
%
% The covariates printed are the raw covariates as entered into SPM, with
% the exception of the global value, which is printed after any grand mean
% scaling.
%
%                           ----------------
%
% FORMAT spm_DesRep('DesMtx',xX,fnames,xs)
% Produces a one-page graphical summary of the design matrix
% xX      - Structure containing design matrix information
%         - the first of {xX.nX, xKXs.X, xX.X} is used for display
% .nX     - Desgin matrix already scaled for display
% .xKXs.X - temporally filtered design matrix (within space structure)
% .X      - "raw" design matrix (as setup by spm{,_fmri}_ui routines)
%
% .Xnames - [optional] px1 CellStr of parameter names
%
% P       - [optional] nx1 CellStr of filenames (i.e. {V.fname}')
% xs      - [optional] structure of extra strings containing descriptive
%          information which is printed at the foot of the page.
%          The field names are used as sub-headings, the field values
%          (which must be strings or CellStr) printed alongside.
%
% The design matrix is labelled with the corresponding parameter and
% file names, and is displayed as an image scaled (using
% spm_DesMtx('sca',...) such that zero is mid-grey, -1 is black, and +1
% is white. Covariates exceeding this randge are scaled to fit.
%
% The design matrix is "clickable": Clicking (and dragging) around the
% design matrix image reports the corresponding value of the design
% matrix ('normal' click - "left" mouse button usually), the image
% filename ('extend' mouse click - "middle" mouse), or parameter name
% ('alt' click - "right" mouse). (spm_DesRep('SurfDesMtx_CB') handles
% the CallBacks)
%
% Under the design matrix the parameter estimability is displayed as a
% 1xp matrix of grey and white squares. Parameters that are not
% uniquely specified by the model are shown with a grey patch. CLicking
% on (and dragging around) the estimability image reports the parameter
% names and their estimability. (spm_DesRep('SurfEstIm_CB') handles
% the CallBacks)
% 
%                           ----------------
%
% FORMAT spm_DesRep('Covs',xC,X,Xnames)
% Plots the covariates and describes how thay are inserted into the model.
% xC     - Covariate structure array (see spm_spm_ui.m for details)
%          ('rcname','rc','descrip','cname' & 'cols' fields used)
% X      - nxp Design matrix
% Xnames - px1 CellStr of parameter names
%
% Covariates are plotted, one per page, overlaid on the design matrix.
% The description strings in the xC covariate structure array are
% displayed. The corresponding design matrix column(s) is(are)
% highlighted.
%
% ======================================================================
% Utility functions and CallBack handlers:
%
% FORMAT s = spm_DesRep('ScanTick',nScan,lim)
% Pares down 1:nScan to at most lim items, showing every 2nd/3rd/... as
% necessary to pair  down to <lim items. Always ends with nScan so
% #images is indicated.
% nScan  - number of scans
% lim    - limit to number of elements of s
% s      - 1:nScan pared down accordingly
%
%                           ----------------
%
% FORMAT spm_DesRep('SurfDesMtx_CB')
% 'ButtonDownFcn' CallBack for surfing clickable design matrix images
% FORMAT spm_DesRep('SurfDesMtxMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing clickable design matrix images
%
% Setting spm_DesRep('SurfDesMtx_CB') as the 'ButtonDownFcn' on a
% design matrix image (and saving design info in it's UserData) allows
% the user to "surf" the design matrix.
%
% "Surfing" (i.e. clicking on the image and dragging around it) with
% the 'normal' mouse button (usually the "left" button) reports desgin
% matrix values, the 'extend' ("middle" button) reports the
% corresponding image file name, and the 'alt' ("right" button" reports
% the parameter names.
%
% The design matrix, parameter names and image file names should be
% saved in the UserData of the image object as a structure with fields
% 'X','Xnames' & 'fnames' respectively. The code is robust to any of
% these being empty or mis-specified - surfing simply reports "no
% cached data".
%
%                           ----------------
%
% FORMAT spm_DesRep('SurfEstIm_CB')
% 'ButtonDownFcn' CallBack for surfing clickable parameter estimability images
% FORMAT spm_DesRep('surfEstImMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing parameter estimability images
%
% Setting spm_DesRep('SurfEstIm_CB') as the 'ButtonDownFcn' on a
% parameter estimability image (and saving design info in it's UserData) allows
% the user to "surf" the parameter estimability image, with the
% parameter at the current point and it's estimability being reported
% continuously.
%
% The binary parameter estimability matrix and parameter names should
% be saved in the UserData of the image object as a structure with
% fields 'est' & 'Xnames' respectively. The code is robust to any of
% these being empty or mis-specified - surfing simply reports "no
% cached data".
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%


%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...'), end

switch lower(varargin{1})

case 'files&factors'                         %-Summarise files & factors
%=======================================================================
% spm_DesRep('Files&Factors',P,I,xC,sF,xs)
if nargin<4, error('insufficient arguments'), end
P  = varargin{2};
I  = varargin{3};
xC = varargin{4};
sF = varargin{5};
if nargin<6, xs=[]; else, xs = varargin{6}; end
			%-Structure of description strings

[P,CPath] = spm_str_manip(P,'c');	%-extract common path component
nScan     = size(I,1);			%-#images
bL        = any(diff(I,1),1); 		%-Multiple factor levels?

%-Get graphics window & window scaling
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph)
FS     = spm('FontSizes');

%-Display header information
%-----------------------------------------------------------------------
hTax = axes('Position',[0.03,0.85,0.94,0.1],...
	'DefaultTextFontSize',FS(9),...
	'XLim',[0,1],'YLim',[0,1],...
	'Visible','off');

text(0.5,1,'Statistical analysis: Image files & covariates...',...
	'Fontsize',FS(14),'Fontweight','Bold',...
	'HorizontalAlignment','center')

dx1 = 0.05;
dx2 = 0.08;

x = 0; text(x+.02,.1,'image #','Rotation',90)
if bL(4), x=x+dx1; text(x+.01,.1,sF{4},'Rotation',90), end
if bL(3), x=x+dx1; text(x+.01,.1,sF{3},'Rotation',90), end
if bL(2), x=x+dx1; text(x+.01,.1,sF{2},'Rotation',90), end
if bL(1), x=x+dx1; text(x+.01,.1,sF{1},'Rotation',90), end

for j = 1:length(xC)
	n = size(xC(j).rc,2);
	if n>1, tmp=xC(j).cname; else, tmp={xC(j).rcname}; end
	for k=1:n
		x=x+dx2;
		text(x,.1,tmp{k},'Rotation',90,'Interpreter','TeX')
	end
end

x=x+dx2;
text(x,0.65,'Base directory:','FontWeight','Bold')
text(x,0.5,CPath,'FontSize',FS(8))
text(x,0.2,'filename tails...')

line('XData',[0 1],'YData',[0 0],'LineWidth',3,'Color','r')

%-Tabulate file & covariate information
%-----------------------------------------------------------------------
hAx = axes('Position',[0.03,0.05,0.94,0.8],...
	'DefaultTextFontSize',FS(8),...
	'Units','points',...
	'Visible','off');
AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])

dy = FS(9); y0 = floor(AxPos(4)) -dy; y  = y0;

for i = 1:nScan

	%-Scan indices
	x = 0; text(x,y,sprintf('%03d',i))
	if bL(4), x=x+dx1; text(x,y,sprintf('%02d',I(i,4))), end
	if bL(3), x=x+dx1; text(x,y,sprintf('%02d',I(i,3))), end
	if bL(2), x=x+dx1; text(x,y,sprintf('%02d',I(i,2))), end
	if bL(1), x=x+dx1; text(x,y,sprintf('%02d',I(i,1))), end

	%-Covariates
	for j = 1:length(xC)
		for k=1:size(xC(j).rc,2)
			x=x+dx2;
			text(x,y,sprintf('%6g',xC(j).rc(i,k)),...
				'HorizontalAlignment','Center')
		end
	end

	%-Filename tail
	x=x+dx2; text(x,y,P{i})

	y=y-dy;

	%-Paginate if necessary
	if y<dy
		text(0.5,0,sprintf('Page %d',spm_figure('#page')),...
			'FontSize',FS(8),'FontAngle','italic')
		spm_figure('NewPage',[hAx;get(hAx,'Children')])
		hAx = axes('Units','points','Position',AxPos,...
			'DefaultTextFontSize',FS(8),'YLim',[0,AxPos(4)],...
			'Visible','off');
		y = y0;
		text(y,0,'continued...','FontAngle','Italic')
	end
end

line('XData',[0 1],'YData',[y y],'LineWidth',3,'Color','r')


%-Display description strings
% (At bottom of current page - hope there's enough room!)
%-----------------------------------------------------------------------
if ~isempty(xs)
	y = y - 2*dy;
	for sf = fieldnames(xs)'
		text(0.3,y,[strrep(sf{1},'_',' '),' :'],...
			'HorizontalAlignment','Right','FontWeight','Bold',...
			'FontSize',FS(9))
		s = getfield(xs,sf{1});
		if ~iscellstr(s), s={s}; end
		for i=1:prod(size(s))
			text(0.31,y,s{i},'FontSize',FS(9))
			y=y-dy;
		end
	end
end

%-Register last page if paginated
if spm_figure('#page')>1
	text(0.5,0,sprintf('Page %d/%d',spm_figure('#page')*[1,1]),...
		'FontSize',FS(8),'FontAngle','italic')
	spm_figure('NewPage',[hAx;get(hAx,'Children')])
end



case 'desmtx'                                    %-Display design matrix
%=======================================================================
% spm_DesRep('DesMtx',xX,fnames,xs)
if nargin<2, error('insufficient arguments'), end
if ~isstruct(varargin{2}), error('design matrix structure required'), end
if nargin<3, fnames={}; else, fnames=varargin{3}; end
if nargin<4, xs=[]; else, xs=varargin{4}; end


%-Locate DesMtx (X), scaled DesMtx (nX) & get parameter names (Xnames)
%-----------------------------------------------------------------------
if isfield(varargin{2},'xKXs') & ...
		~isempty(varargin{2}.xKXs) & isstruct(varargin{2}.xKXs)
	iX = 1;
	[nScan,nPar] = size(varargin{2}.xKXs.X);
elseif isfield(varargin{2},'X') & ~isempty(varargin{2}.X)
	iX = 0;
	[nScan,nPar] = size(varargin{2}.X);
else
	error('Can''t find DesMtx in this structure!')
end

if isfield(varargin{2},'nKX') & ~isempty(varargin{2}.nKX)
	inX = 1; else, inX = 0; end

if isfield(varargin{2},'Xnames') & ~isempty(varargin{2}.Xnames)
	Xnames = varargin{2}.Xnames; else, Xnames = {}; end


%-Get graphics window & FontSizes
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)
FS = spm('FontSizes');

%-Title
%-----------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
	'DefaultTextFontSize',FS(9),...
	'XLim',[0,1],'YLim',[0,1],...
	'Visible','off');

text(0.5,0.95,'Statistical analysis: Design',...
	'Fontsize',FS(14),'Fontweight','Bold',...
	'HorizontalAlignment','center')

line('XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r')
line('XData',[0.3 0.7],'YData',[0.28 0.28],'LineWidth',3,'Color','r')


%-Display
%-----------------------------------------------------------------------
hDesMtx   = axes('Position',[.07 .4 .6 .4]);
if inX		%-Got a scaled DesMtx
	hDesMtxIm = image((varargin{2}.nKX + 1)*32);
elseif iX	%-No scaled DesMtx, DesMtx in .xKXs structure
	hDesMtxIm = image((spm_DesMtx('sca',varargin{2}.xKXs.X,Xnames) + 1)*32);
else		%-No scaled DesMtx, no .xKXs, DesMtx in .X
	hDesMtxIm = image((spm_DesMtx('sca',varargin{2}.X,     Xnames) + 1)*32);
end
set(hDesMtx,'TickDir','out')
ylabel('images')
xlabel('parameters')

%-Parameter names
if ~isempty(Xnames)
	axes('Position',[.07 .8 .6 .1],'Visible','off',...
		'DefaultTextFontSize',FS(8),'DefaultTextInterpreter','TeX',...
		'XLim',[0,nPar]+0.5)
	for i = 1:nPar, text(i,.05,Xnames{i},'Rotation',90), end
end

%-Filenames
% ( Show at most 32, showing every 2nd/3rd/4th/... as necessary to pair )
% ( down to <32 items. Always show last item so #images is indicated.   )     
if ~isempty(fnames)
	s = spm_DesRep('ScanTick',nScan,32);
	set(hDesMtx,'YTick',s)
	axes('Position',[.68 .4 .3 .4],'Visible','off',...
		'DefaultTextFontSize',FS(8),...
		'YLim',[0,nScan]+0.5,'YDir','Reverse')
	for i=s, text(0,i,spm_str_manip(fnames{i},'a40')), end
end

%-Setup callbacks to allow interrogation of design matrix
%-----------------------------------------------------------------------
if iX, 	set(hDesMtxIm,'UserData',...
	struct('X',varargin{2}.xKXs.X,'Xnames',{Xnames},'fnames',{fnames}))
else, 	set(hDesMtxIm,'UserData',...
	struct('X',varargin{2}.X,     'Xnames',{Xnames},'fnames',{fnames}))
end
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')

%-Parameter estimability/uniqueness
%-----------------------------------------------------------------------
tmp = get(Fgraph,'DefaultTextInterpreter');
set(Fgraph,'DefaultTextInterpreter','TeX');
hPEstAx   = axes('Position',[.07 .325 .6 .025],'DefaultTextInterpreter','TeX');
if iX,	est = spm_SpUtil('IsCon',varargin{2}.xKXs);
else,	est = spm_SpUtil('IsCon',varargin{2}.X); end
hParEstIm = image((est+1)*32);
set(hPEstAx,...
	'XLim',[0,nPar]+.5,'XTick',[1:nPar-1]+.5,'XTickLabel','',...
	'YLim',[0,1]+.5,'YDir','reverse','YTick',[],...
	'Box','on','TickDir','in','XGrid','on','GridLineStyle','-');
xlabel('paremeter estimability')
text((nPar+0.5 + nPar/30),1,...
	'(gray \rightarrow \beta not uniquely specified)',...
	'Interpreter','TeX','FontSize',FS(8))
set(hParEstIm,'UserData',struct('est',est,'Xnames',{Xnames}))
set(hParEstIm,'ButtonDownFcn','spm_DesRep(''SurfEstIm_CB'')')
set(Fgraph,'DefaultTextInterpreter',tmp);


%-Design descriptions
%-----------------------------------------------------------------------
if ~isempty(xs)
	hAx = axes('Position',[0.03,0.05,0.94,0.22],'Visible','off');
	
	set(hAx,'Units','points');
	AxPos = get(hAx,'Position');
	set(hAx,'YLim',[0,AxPos(4)])
	
	dy = FS(9); y0 = floor(AxPos(4)) -dy; y = y0;

	text(0.3,y,'Design description...',...
		'HorizontalAlignment','Center',...
		'FontWeight','Bold','FontSize',FS(11))
	y=y-2*dy;
	
	for sf = fieldnames(xs)'
		text(0.3,y,[strrep(sf{1},'_',' '),' :'],...
			'HorizontalAlignment','Right','FontWeight','Bold',...
			'FontSize',FS(9))
		s = getfield(xs,sf{1});
		if ~iscellstr(s), s={s}; end
		for i=1:prod(size(s))
			text(0.31,y,s{i},'FontSize',FS(9))
			y=y-dy;
		end
	end
end



case 'covs'                %-Plot and describe covariates (one per page)
%=======================================================================
% spm_DesRep('Covs',xC,X,Xnames)
if nargin<4, error('insufficient arguments'), end
xC     = varargin{2};	%-Struct array of covariate information
X      = varargin{3};	%-nxp Design matrix
Xnames = varargin{4};	%-px1 CellStr of parameter names

if ~length(xC), return, end

%-Get graphics window & window scaling
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph)
FS = spm('FontSizes');

%-Title
%-----------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
	'DefaultTextFontSize',FS(9),...
	'XLim',[0,1],'YLim',[0,1],...
	'Visible','off');

text(0.5,0.95,'Statistical analysis: Covariates',...
	'Fontsize',FS(14),'Fontweight','Bold',...
	'HorizontalAlignment','center')

text(0.5,0.82,'(covariates plotted over transposed design matrix)',...
	'FontSize',FS(8),'HorizontalAlignment','center')

line('XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r')
line('XData',[0.3 0.7],'YData',[0.44 0.44],'LineWidth',3,'Color','r')


%-Design matrix (as underlay for plots) and parameter names
%-----------------------------------------------------------------------
%-Design matrix
nX = spm_DesMtx('sca',X,Xnames);
hDesMtx = axes('Position',[.1 .5 .7 .3]);
image((nX' + 1)*32);
set(hDesMtx,'Visible','off')

%-Parameter names
hParAx = axes('Position',[.8 .5 .2 .3],'Visible','off',...
	'DefaultTextFontSize',FS(8),'DefaultTextInterpreter','TeX',...
	'YLim',[0.5,size(X,2)+0.5],'YDir','Reverse');
hPNames = zeros(size(X,2),1);
for i = 1:size(X,2), hPNames(i) = text(.05,i,Xnames{i}); end


%-Covariates - one page each
%-----------------------------------------------------------------------
nScan = size(X,1);
for i = 1:length(xC)

	%-Title
	%---------------------------------------------------------------
	hSTitle = text(0.5,0.87,sprintf('%d : %s',i,xC(i).rcname),...
			'Parent',hTax,...
			'HorizontalAlignment','center',...
			'FontSize',FS(13),'FontWeight','Bold');

	%-Plot
	%---------------------------------------------------------------
	hAx = axes('Position',[.1 .5 .7 .3],...
			'TickDir','out','Box','off','Color','none',...
			'NextPlot','add',...
			'XLim',[0,nScan]+0.5);
	plot(xC(i).rc,'LineWidth',2)
	if nScan<48, plot(xC(i).rc,'.k','MarkerSize',20); end
	xlabel('image #')
	ylabel('covariate value')


	%-Descriptions
	%---------------------------------------------------------------
	hDAx = axes('Position',[0.03,0.1,0.94,0.30],'Visible','off');
	
	set(hDAx,'Units','points');
	tmp = get(hDAx,'Position');
	set(hDAx,'YLim',[0,tmp(4)])
	
	dy = FS(9); y0 = floor(tmp(4)) -dy; y = y0;

	%-Description strings from xC(i).descrip
	text(0.3,y,'Details :',...
		'HorizontalAlignment','Right',...
		'FontWeight','Bold','FontSize',FS(9))
	s = xC(i).descrip;
	if ~iscellstr(s), s={s}; end
	for j=1:prod(size(s))
		text(0.31,y,s{j},'FontSize',FS(9))
		y=y-dy;
	end
	y=y-dy;

	%-Key (if block of covariates entered)
	%---------------------------------------------------------------
	if size(xC(i).rc,2)>1
		ColorOrder = get(hAx,'ColorOrder');
		text(0.3,y,'Key :',...
			'HorizontalAlignment','Right',...
			'FontWeight','Bold','FontSize',FS(9))
		for j = 1:size(xC(i).rc,2)
			color = ColorOrder(mod(j-1,size(ColorOrder,1))+1,:);
			if size(xC(i).rc,2)==length(xC(i).cname)
				str = xC(i).cname{j};
			else
				str = sprintf('column %d',j);
			end
			text(0.31,y,str,'FontSize',FS(9),...
				'Color',color)
			text(0.5,xC(i).rc(1,j),[str,' \rightarrow'],...
				'Parent',hAx,...
				'FontSize',FS(8),'FontWeight','Bold',...
				'HorizontalAlignment','Right',...
				'Interpreter','TeX',...
				'Color',color)
			y=y-dy;
		end
		y=y-dy;
	end


	%-Associated parameters
	%---------------------------------------------------------------
	text(0.3,y,'Design matrix columns :',...
		'HorizontalAlignment','Right',...
		'FontWeight','Bold','FontSize',FS(9))
	if isempty(xC(i).cols)
		text(0.31,y,'(none)','FontSize',FS(9))
	else
		for j = xC(i).cols
			text(0.31,y,sprintf('%d : %s',j,Xnames{j}),...
				'FontSize',FS(9),'Interpreter','TeX')
			y=y-dy;
		end
	end
	y=y-dy;


	%-Highlight parameter names
	%---------------------------------------------------------------
	hCurPNames = hPNames(xC(i).cols);
	set(hCurPNames,'Color','r','FontWeight','Bold','FontSize',FS(8))


	%-Paginate (if more than one covariate)
	%---------------------------------------------------------------
	if length(xC)>1
		spm_figure('NewPage',[hSTitle; hAx; get(hAx,'Children');...
			hCurPNames; hDAx; get(hDAx,'Children')]);
	end

end


case 'scantick'
%=======================================================================
% spm_DesRep('ScanTick',nScan,lim)
% ( Show at most 32, showing every 2nd/3rd/4th/... as necessary to pair )
% ( down to <32 items. Always show last item so #images is indicated.    )     
if nargin<3, lim=32; else, lim=varargin{3}; end
if nargin<2, error('insufficient arguments'), end
nScan = varargin{2};

p = max(1,ceil(nScan/lim));
s = 1:p:nScan; s(end)=nScan;

varargout = {s,lim};


case {'surfdesmtx_cb','surfdesmtxmo_cb'}  %-Callbacks for surfing DesMtx
%=======================================================================
% spm_DesRep('SurfDesMtx_CB')
% spm_DesRep('SurfDesMtxMo_CB')

if strcmp(lower(varargin{1}),'surfdesmtx_cb')
	set(0,'UserData',get(get(gca,'XLabel'),'String'))
	set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfDesMtxMo_CB'')')
	set(gcbf,'WindowButtonUpFcn',[...
		'set(get(gca,''XLabel''),''String'',get(0,''UserData'')),',...
		'set(0,''UserData'',[]),',...
		'set(gcbf,''WindowButtonMotionFcn'',''''),',...
		'set(gcbf,''WindowButtonUpFcn'','''')'])
end

mm  = [get(gca,'YLim')',get(gca,'XLim')']+[.5,.5;-.5,-.5];
ij  = get(gca,'CurrentPoint');
ij  = round(min(max(ij(1,[2,1]),mm(1,:)),mm(2,:)));

switch get(gcbf,'SelectionType')
case 'normal'
	try, str = sprintf('X(%d,%d) = %g',ij(1),ij(2),...
		subsref(get(gco,'UserData'),...
		struct('type',{'.','()'},'subs',{'X',{ij(1),ij(2)}})));
	catch, str='(no cached design matrix to surf)'; end
case 'extend'
	try, str = sprintf('Image %d: %s',ij(1),...
		spm_str_manip(...
		subsref(get(gco,'UserData'),...
		struct('type',{'.','{}'},'subs',{'fnames',{ij(1)}})),'a40'));
	catch, str='(no cached image filenames to surf)'; end
case 'alt'
	try, str = sprintf('Parameter %d: %s',ij(2),...
		subsref(get(gco,'UserData'),...
		struct('type',{'.','{}'},'subs',{'Xnames',{ij(2)}})));
	catch, str='(no cached parameter names to surf)'; end
case 'open'
	return
end

set(get(gca,'XLabel'),'String',str)


case {'surfestim_cb','surfestimmo_cb'}  %-Callbacks for surfing ParEstIm
%=======================================================================
% spm_DesRep('SurfEstIm_CB')
% spm_DesRep('surfEstImMo_CB')

if strcmp(lower(varargin{1}),'surfestim_cb')
	set(0,'UserData',get(get(gca,'XLabel'),'String'))
	set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''surfEstImMo_CB'')')
	set(gcbf,'WindowButtonUpFcn',[...
		'set(get(gca,''XLabel''),''String'',get(0,''UserData'')),',...
		'set(0,''UserData'',[]),',...
		'set(gcbf,''WindowButtonMotionFcn'',''''),',...
		'set(gcbf,''WindowButtonUpFcn'','''')'])
end

mm  = [get(gca,'XLim')]+[.5,-.5];
i   = get(gca,'CurrentPoint');
i   = round(min(max(i(1,1),mm(1)),mm(2)));

switch get(gcbf,'SelectionType')
case 'normal'
	try, tmp = {' (not unique)',' (unique)'};
	str = sprintf('Parameter %d : %s%s',...
		i,...
		subsref(get(gco,'UserData'),...
			struct('type',{'.','{}'},'subs',{'Xnames',{i}})),...
		tmp{subsref(get(gco,'UserData'),...
			struct('type',{'.','()'},'subs',{'est',{i}}))+1});
	catch, str='(no cached data to surf)'; end
case {'extend','alt','open'}
	return
end

set(get(gca,'XLabel'),'String',str)



otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string')



%=======================================================================
end
