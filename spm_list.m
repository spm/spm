function spm_list(SPM,VOL,Dis,Num,title)
% Display and analysis of SPM{.}
% FORMAT spm_list(SPM,VOL,Dis,Num,title)
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
%
% VOL    - structure containing details of volume analysed
%        - required fields are:
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
%
% Dis    - Minimum distance between maxima                 {default = 8mm}
% Num    - Maxiumum number of maxima tabulated per cluster {default = 2}
% title  - title text for table
%
% (see spm_getSPM for further details of SPM & VOL structures)
%_______________________________________________________________________
%
% spm_list characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected subsets
% or clusters).  The significance of the results are based on set, cluster
% and voxel-level inferences using distributional approximations from the
% Theory of Gaussian Feilds.  These distributions assume that the SPM is
% a reasonable lattice approximation with known component field smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C = P(u,k,C).  For each cluster
% of size K the cluster-level inference is based on P(u,K,1) and for each
% voxel (or selected maxima) of height U, in that cluster, the voxel-level
% inference is based on P(U,0,1).  All three levels of inference are
% supported with a tabular presentation of the p values and the underlying
% statistic (u, k or c).  The table is grouped by regions and sorted on
% the Z-variate of the primary maxima.  See 'Sections' in the help facility
% and spm_maxima.m for a more complete characterization of maxima within
% a region. Z-variates (based on the uncorrected p value) are the Z score
% equivalent of the statistic. Volumes are expressed in voxels.
%
% Click on co-ordinates to jump results cursor to location ****
% Click on tabulated data to extract accurate value to MatLab
% workspace. (Look in the MatLab command window.)
%
%_______________________________________________________________________
% %W% Karl Friston %E%


%-Parse arguments
%-----------------------------------------------------------------------
if nargin < 2,   error('insufficient arguments'), end
if nargin < 5,   title = spm_str_manip(SPM.swd,'a50'); end
if nargin < 4,   Num = []; end
if isempty(Num), Num = 3;  end
if nargin < 3,   Dis = []; end
if isempty(Dis), Dis = 8;  end


%-Setup graphics pane
%-----------------------------------------------------------------------
spm('Pointer','Watch')
Fgraph    = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph)
FS        = spm('FontSizes');


%-Characterize excursion set in terms of maxima
% (sorted on Z values and grouped by regions)
%-----------------------------------------------------------------------
n         = SPM.n;
STAT      = SPM.STAT;
df        = SPM.df;
u         = SPM.u;
k         = SPM.k;
R         = VOL.R;

if ~length(SPM.Z)
	spm('Pointer','Arrow')
	msgbox('No voxels survive masking & threshold(s)!',...
		sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'help','modal')
	return
end

[N Z XYZ A] = spm_max(SPM.Z,SPM.XYZ);

%-Convert cluster sizes and extent threshold from voxels to resels
%-----------------------------------------------------------------------
v2r       = 1/prod(VOL.FWHM);				% voxels to resels
N         = N*v2r;
k         = k*v2r;

%-Convert maxima locations from voxels to mm
%-----------------------------------------------------------------------
XYZmm     = VOL.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];


%-Table axes & headings
%=======================================================================
hAx   = axes('Position',[0.05 0.1 0.9 0.4],...
	'DefaultTextFontSize',FS(8),...
	'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline',...
	'Units','points',...
	'Visible','off');
AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
dy    = FS(9);
y     = floor(AxPos(4)) - dy;

text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',title],...
	'FontSize',FS(11),'FontWeight','Bold');	y = y - dy/2;
line([0 1],[y y],'LineWidth',3,'Color','r'),	y = y - 9*dy/8;

%-Construct table header
%-----------------------------------------------------------------------
set(gca,'DefaultTextFontName','Helvetica','DefaultTextFontSize',FS(8))

if max(A) > 1

	text(0.01,y,		'set-level','FontSize',FS(9))
	line([0.00,0.11],[y-dy/4,y-dy/4],'LineWidth',0.5,'Color','r')
	text(0.02,y-9*dy/8,	'\itp ')
	text(0.09,y-9*dy/8,	'\itc ')
end

text(0.22,y,		'cluster-level','FontSize',FS(9))
line([0.16,0.42],[y-dy/4,y-dy/4],'LineWidth',0.5,'Color','r')
text(0.16,y-9*dy/8,	'\itp \rm_{corrected}')
text(0.28,y-9*dy/8,	'\itk')
text(0.33,y-9*dy/8,	'\itp \rm_{uncorrected}')

text(0.60,y,		'voxel-level','FontSize',FS(9))
line([0.48,0.83],[y-dy/4,y-dy/4],'LineWidth',0.5,'Color','r')
text(0.48,y-9*dy/8,	'\itp \rm_{corrected}')
text(0.60,y-9*dy/8,	sprintf('\\it%c',STAT))
text(0.68,y-9*dy/8,	'(\itZ\rm_\equiv)')
text(0.75,y-9*dy/8,	'\itp \rm_{uncorrected}')

text(0.90,y - dy/2,['x,y,z \fontsize{',num2str(FS(8)),'}\{mm\}']);

%-Headers for text table...
TabDat.tit = title;
TabDat.hdr = {	'set',		'c';...
		'set',		'p';...
		'cluster',	'p (cor)';...
		'cluster',	'k';...
		'cluster',	'p (unc)';...
		'voxel',	'p (cor)';...
		'voxel',	STAT;...
		'voxel',	'equivZ';...
		'voxel',	'p (unc)';...
		'',		'x,y,z {mm}'}';

y     = y - 7*dy/4;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - 5*dy/4;

%-Pagination variables
%-----------------------------------------------------------------------
y0    = y;
hPage = [];

set(gca,'DefaultTextFontName','Courier','DefaultTextFontSize',FS(7))

%-Table proper (& note all data in cell array)
%=======================================================================

%-Set-level p values {c}
%-----------------------------------------------------------------------

c     = max(A);					%-Number of clusters
Pc    = spm_P(c,k,u,df,STAT,R,n);		%-Set-level p-value

% do not display if reporting a single cluster
%-----------------------------------------------------------------------
if c > 1;

	h     = text(0.00,y,sprintf('%-0.3f',Pc),'FontWeight','Bold',...
		'UserData',Pc,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.08,y,sprintf('%g',c),'FontWeight','Bold',...
		'UserData',c,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

end
TabDat = {Pc,c};				%-Table data
TabLin = 1;					%-Table data line

%-Local maxima p-values & statistics
%=======================================================================
while max(Z)

	% Paginate if necessary
	%---------------------------------------------------------------
	if y < (Num + 1)*dy
		h     = text(0.5,-5*dy,...
			sprintf('Page %d',spm_figure('#page')),...
			'FontName','Helvetica','FontAngle','Italic',...
			'FontSize',FS(8));
		spm_figure('NewPage',[hPage,h])
		hPage = [];
		y     = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U,i]   = max(Z);			% largest maxima
	j       = find(A == A(i));		% maxima in cluster


    	%-Compute cluster {k} and voxel-level {u} p values for this cluster
    	%---------------------------------------------------------------
	Pz      = spm_P(1,0,   U,df,STAT,1,n);	% uncorrected p value
	Pu      = spm_P(1,0,   U,df,STAT,R,n);	% corrected     {based on Z)
	[Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n);	% [un]corrected {based on k)
	Nv      = N(i)/v2r;			% extent        {voxels}
	Ze      = spm_invNcdf(1 -Pz);		% Equivalent Z-variate


	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
	h     = text(0.17,y,sprintf('%0.3f',Pk),	'FontWeight','Bold',...
		'UserData',Pk,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.27,y,sprintf('%0.0f',Nv),	'FontWeight','Bold',...
		'UserData',N(i),'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.35,y,sprintf('%0.3f',Pn),	'FontWeight','Bold',...
		'UserData',Pn,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

	h     = text(0.49,y,sprintf('%0.3f',Pu),	'FontWeight','Bold',...
		'UserData',Pu,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.57,y,sprintf('%6.2f',U),		'FontWeight','Bold',...
		'UserData',U,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.66,y,sprintf('(%5.2f)',Ze),	'FontWeight','Bold',...
		'UserData',Ze,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(0.77,y,sprintf('%0.3f',Pz),	'FontWeight','Bold',...
		'UserData',Pz,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

	h     = text(0.88,y,sprintf('%3.0f %3.0f %3.0f',XYZmm(:,i)),...
		'FontWeight','Bold',...
		'ButtonDownFcn',...
		'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''));',...
		'Interruptible','off','BusyAction','Cancel',...
		'UserData',XYZmm(:,i));
	hPage = [hPage, h];
 
	y     = y - dy;
	
	[TabDat{TabLin,3:10}] = deal(Pk,Nv,Pn,Pu,U,Ze,Pz,XYZmm(:,i));
	TabLin = TabLin + 1;

	%-Print Num secondary maxima (> Dis mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-Z(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
	    d = j(q(i));
	    if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;
		if length(D) < Num
			
			% voxel-level p values {Z}
			%-----------------------------------------------
			Pz    = spm_P(1,0,Z(d),df,STAT,1,n);
			Pu    = spm_P(1,0,Z(d),df,STAT,R,n);
			Ze    = spm_invNcdf(1 - Pz);

			h     = text(0.49,y,sprintf('%0.3f',Pu),...
				'UserData',Pu,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(0.57,y,sprintf('%6.2f',Z(d)),...
				'UserData',Z(d),...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(0.66,y,sprintf('(%5.2f)',Ze),...
				'UserData',Ze,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(0.77,y,sprintf('%0.3f',Pz),...
				'UserData',Pz,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(0.88,y,...
				sprintf('%3.0f %3.0f %3.0f',XYZmm(:,d)),...
				'ButtonDownFcn',[...
					'spm_mip_ui(''SetCoords'',',...
					'get(gcbo,''UserData''));'],...
				'Interruptible','off','BusyAction','Cancel',...
				'UserData',XYZmm(:,d));
			hPage = [hPage, h];
			D     = [D d];
			y     = y - dy;
			[TabDat{TabLin,6:10}] = deal(Pu,Z(d),Ze,Pz,XYZmm(:,d));
			TabLin = TabLin+1;
		end
	    end
	end
	Z(j) = Z(j)*0;		% Zero local maxima for this cluster
end				% end region


%-Number and register last page (if paginated)
%-----------------------------------------------------------------------
if spm_figure('#page')>1
	h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page')*[1,1]),...
		'FontName','Helvetica','FontSize',FS(8),'FontAngle','Italic');
	spm_figure('NewPage',[hPage,h])
end

%-Table filtering note
%===========================================================================
str = sprintf(['table shows at most %d subsidiary maxima ',...
	'> %.1fmm apart per cluster'],Num,Dis);
text(0.5,4,str,'HorizontalAlignment','Center','FontName','Helvetica',...
    'FontSize',FS(8),'FontAngle','Italic')


%-Volume, resels and smoothness 
%===========================================================================
FWHMmm          = VOL.FWHM.*VOL.VOX'; 				% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n);
Pu              = spm_P(1,0,u,df,STAT,R,n);
[P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n);

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
line([0 1],[0 0],'LineWidth',1,'Color','r')
set(gca,'DefaultTextFontName','Helvetica',...
	'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
TabFut = cell(4,2);
TabFut{1} = ...
	sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f corrected)',...
		 STAT,u,Pz,Pu);
TabFut{2} = ...
sprintf('Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f corrected)',...
	         k/v2r,Pn,P);
TabFut{3} = ...
	sprintf('Expected voxels per cluster, <k> = %0.3f',En/v2r);
TabFut{4} = ...
	sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
TabFut{5} = ...
	sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
TabFut{6} = ...
	sprintf(['Smoothness FWHM = %0.1f %0.1f %0.1f {mm} ',...
		 ' = %0.1f %0.1f %0.1f {voxels}'],FWHMmm,VOL.FWHM);
TabFut{7} = ...
	sprintf(['Volume: S = %0.0f voxels = %0.2f resels ',...
	         '(1 resel = %0.1f voxels)'],VOL.S,R(end),prod(VOL.FWHM));
TabFut{8} = ...
	sprintf('');

text(0.0,-1*dy,TabFut{1},...
	'UserData',[u,Pz,Pu],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-2*dy,TabFut{2},...
	'UserData',[k/v2r,Pn,P],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-3*dy,TabFut{3},...
	'UserData',En/v2r,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-4*dy,TabFut{4},...
	'UserData',Em*Pn,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-1*dy,TabFut{5},...
	'UserData',df,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-2*dy,TabFut{6},...
	'UserData',FWHMmm,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-3*dy,TabFut{7},...
	'UserData',[VOL.S,R(end),prod(VOL.FWHM)],...
	'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-4*dy,TabFut{8},...
	'UserData',[],'ButtonDownFcn','get(gcbo,''UserData'')')


%-Store TabDat in UserData of axes
%=======================================================================
% TabDat, TabHdr, TabFut
h      = uicontextmenu('Tag','TabDat');
set(gca,'UIContextMenu',h,...
	'Visible','on',...
	'XColor','w','YColor','w')
uimenu(h,'Label','Extract table data in cell array',...
	'Tag','TD_Xdat',...
	'CallBack','',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','on','Label','Print table as text',...
	'Tag','TD_TxtTab',...
	'CallBack','',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','on','Label','Print table (to cluster level) as text',...
	'Tag','TD_TxtTab',...
	'CallBack','',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','on','Label','Print table (to voxel level) as text',...
	'Tag','TD_TxtTab',...
	'CallBack','',...
	'Interruptible','off','BusyAction','Cancel');



%-End
%=======================================================================
spm('Pointer','Arrow')
