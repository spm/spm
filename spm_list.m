function spm_list(SPM,VOL,Dis,Num,title)
% Display and analysis of SPM{.}
% FORMAT spm_list(SPM,VOL)
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
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
% equivalent of the statistic. Volumes are expressed in resels.
%
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Default arguments
%-----------------------------------------------------------------------
if nargin<2, error('insufficient arguments'), end
if nargin<5, title=spm_str_manip(SPM.swd,'a50'); end
if nargin<4, Num=[]; end
if isempty(Num), Num=3; end
if nargin<3, Dis=[]; end
if isempty(Dis), Dis=8; end


%-Setup graphics pane
%-----------------------------------------------------------------------
spm('Pointer','Watch')
Fgraph    = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph)
FS        = spm_figure('FontSizes');


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

%-Convert cluster sizes from voxels to resels
N         = N/prod(VOL.FWHM);
%-Convert maxima locations from voxels to mm
XYZmm     = VOL.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];


%-Table axes & headings
%=======================================================================
hAx   = axes('Position',[0.05 0.1 0.9 0.42],...
	'DefaultTextFontSize',FS(1),...
	'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline',...
	'Units','points',...
	'Visible','off');
AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
dy    = FS(2);
y     = floor(AxPos(4)) - dy;

text(0,y,['Statistics:  \it\fontsize{',num2str(FS(2)),'}',title],...
	'FontSize',FS(3),'FontWeight','Bold');	y = y - dy/2;
line([0 1],[y y],'LineWidth',3,'Color','r'),	y = y - 3*dy/2;

%-Construct tables
%-----------------------------------------------------------------------
text(0.00,y,'set-level \{\it{c}\rm\}','FontSize',FS(2));
text(0.15,y,'cluster-level \{\it{k}\rm_{max}\}','FontSize',FS(2));
text(0.38,y,['voxel-level \{\it{',STAT,'}\rm_{max} \equiv \it{Z}\rm_{max}\}'],...
							'FontSize',FS(2));
text(0.66,y,['uncorrected \it{k}\rm & \it{',STAT,'}'],'FontSize',FS(2));
text(0.90,y,['x,y,z \fontsize{',num2str(FS(1)),'}\{mm\}'],'FontSize',FS(2));

y     = y - dy/2;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - 3*dy/2;

%-Pagination variables
%-----------------------------------------------------------------------
y0    = y;
hPage = [];

set(gca,'DefaultTextFontName','Courier','DefaultTextFontSize',FS(1)-1)

%-Set-level p values {c}
%-----------------------------------------------------------------------
c     = max(A);					%-Number of clusters
Pc    = spm_P(c,k,u,df,STAT,R,n);		%-Set-level p-value
str   = sprintf('%-0.3f  (%i)',Pc,c);
h     = text(0.00,y,str,'FontWeight','Bold');
hPage = [hPage, h];


%-Local maxima p-values & statistics
%=======================================================================
while max(Z)

	% Paginate if necessary
	%---------------------------------------------------------------
	if y < (Num+1)*dy
		h     = text(0.5,-5*dy,...
			sprintf('Page %d',spm_figure('#page')),...
			'FontName','Helvetica','FontAngle','Italic',...
			'FontSize',FS(1));
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
	Ze      = spm_invNcdf(1 - Pz);		% Equivalent Z-variate


	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
	h     = text(0.17,y,sprintf('%0.3f  (%0.2f)',Pk,N(i)),...
			'FontWeight','Bold');
	hPage = [hPage, h];
	h     = text(0.38,y,sprintf('%5.3f   (%6.2f \\equiv %5.2f)',...
			Pu,U,Ze),...
			'FontWeight','Bold');
	hPage = [hPage, h];
	h     = text(0.68,y,sprintf('%0.3f',Pn),'FontWeight','Bold');
	hPage = [hPage, h];
	h     = text(0.76,y,sprintf('%0.3f',Pz),'FontWeight','Bold');
	hPage = [hPage, h];
	h     = text(0.88,y,sprintf('%3.0f %3.0f %3.0f',XYZmm(:,i)),...
		'FontWeight','Bold',...
		'ButtonDownFcn','spm_mip_ui(''ShowGreens'')',...
		'Interruptible','off',...
		'UserData',XYZmm(:,i));
	hPage = [hPage, h];
 
	y     = y - dy;

	%-Print Num secondary maxima (> Dis mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-Z(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
	    d    = j(q(i));
	    if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;
		if length(D) < Num
			
			% voxel-level p values {Z}
			%-----------------------------------------------
			Pz    = spm_P(1,0,Z(d),df,STAT,1,n);
			Pu    = spm_P(1,0,Z(d),df,STAT,R,n);
			Ze    = spm_invNcdf(1 - Pz);

			h     = text(0.38,y,...
				sprintf('%5.3f   (%6.2f \\equiv %5.2f)',...
					Pu,Z(d),Ze));
			hPage = [hPage, h];
			h     = text(0.76,y,sprintf('%0.3f',Pz));
			hPage = [hPage, h];
			h     = text(0.88,y,...
				sprintf('%3.0f %3.0f %3.0f',XYZmm(:,d)),...
				'ButtonDownFcn',...
				'spm_mip_ui(''ShowGreens'')',...
				'Interruptible','off',...
				'UserData',XYZmm(:,d));
			hPage = [hPage, h];
			D     = [D d];
			y     = y - dy;
		end
	    end
	end
	Z(j) = Z(j)*0;		% Zero local maxima for this cluster
end				% end region


%-Number and register last page (if paginated)
%-----------------------------------------------------------------------
if spm_figure('#page')>1
	h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page')*[1,1]),...
		'FontName','Helvetica','FontSize',FS(1),'FontAngle','Italic');
	spm_figure('NewPage',[hPage,h])
end


%-Table filtering note
%-----------------------------------------------------------------------
str = sprintf(['table shows at most %d subsidiary maxima ',...
	'>%.1fmm apart per cluster'],Num,Dis);
text(0.5,4,str,'HorizontalAlignment','Center',...
	'FontAngle','Italic','FontSize',FS(1)-1)


%-Volume, resels and smoothness 
%===========================================================================
FWHMmm          = VOL.FWHM.*VOL.VOX'; 				% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n);
[P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n);

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
line([0 1],[0 0],'LineWidth',1,'Color','r')
set(gca,'DefaultTextFontName','Helvetica',...
	'DefaultTextInterpreter','None','DefaultTextFontSize',FS(1))
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f (%0.3f)',u,Pz,P);
text(0.0,-1*dy,str);
str = sprintf('Extent threshold {k} = %0.2f resels, p = %0.3f',k,Pn);
text(0.0,-2*dy,str);
str = sprintf('Expected resels per cluster, E{n} = %0.3f',En);
text(0.0,-3*dy,str);
str = sprintf('Expected number of clusters, E{m} = %0.2f',Em*Pn);
text(0.0,-4*dy,str);

str = sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
text(0.5,-1*dy,str);
str = sprintf(['Smoothness FWHM = %0.1f %0.1f %0.1f {mm} ',...
		' = %0.1f %0.1f %0.1f {voxels}'],FWHMmm,VOL.FWHM);
text(0.5,-2*dy,str);
str = sprintf(['Volume {S} = %0.0f voxels = %0.2f resels ',...
		'(1 resel = %0.1f voxels)'],VOL.S,R(end),prod(VOL.FWHM));
text(0.5,-3*dy,str);
str = sprintf('');
text(0.5,-4*dy,str);

%-End
%=======================================================================
spm('Pointer','Arrow')
