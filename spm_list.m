function spm_list(SPM,VOL,Dis,Num)
% display and analysis of SPM{Z} [SPM{t}]
% FORMAT spm_list(SPM,VOL)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% Dis  - Minimum distance between maxima                 {default = 8mm}
% Num  - Maxiumum number of maxima tabulated per cluster {default = 2}
%
% see spm_getSPM for details
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
global CWD
if nargin<4, Num=3; end
if nargin<3, Dis=8; end

%-Setup graphics pane
%-----------------------------------------------------------------------
spm('Pointer','Watch')
Fgraph    = spm_figure('GetWin','Graphics');
spm_results_ui('ClearPane',Fgraph)
figure(Fgraph)
FS        = spm_figure('FontSizes');


%-Characterize excursion set in terms of maxima
% (sorted on Z values and grouped by regions)
%-----------------------------------------------------------------------
Z         = SPM.Z;
n         = SPM.n;
STAT      = SPM.STAT;
df        = SPM.df;
u         = SPM.u;
k         = SPM.k;
R         = VOL.R;
FWHM      = VOL.FWHM;
[N T M A] = spm_max(Z,VOL.XYZ,VOL.VOX);
N         = N/prod(FWHM); % voxels -> resels


%-Table axes & headings
%=======================================================================
hAx   = axes('Position',[0.1 0.1 0.8 0.42],...
	'DefaultTextFontSize',FS(1),...
	'Units','points',...
	'Visible','off');
AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
dy    = FS(2);
y     = floor(AxPos(4)) - dy;

text(0,y,['Statistics:  ',spm_str_manip(CWD,'a50')],...
	'FontSize',FS(3),'FontWeight','Bold');	y = y - dy;
line([0 1],[y y],'LineWidth',3,'Color','r'),	y = y - dy;

%-Construct tables
%-----------------------------------------------------------------------
text(0.00,y,'set-level {c}'      ,'FontSize',FS(2));
text(0.18,y,'cluster-level {k}  ','FontSize',FS(2));
text(0.42,y,'voxel-level {Z}'    ,'FontSize',FS(2));
text(0.62,y,'uncorrected k & Z'  ,'FontSize',FS(2));
text(0.86,y,'x,y,z {mm}'         ,'FontSize',FS(2));

y     = y - dy;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - dy;

%-Pagination variables
%-----------------------------------------------------------------------
y0    = y;
hPage = [];


%-Set-level p values {c}
%-----------------------------------------------------------------------
c     = max(A);					% number of clusters
Pc    = spm_P(c,k,u,df,STAT,R,n);		% set-level p value

str   = sprintf('%-0.3f   (%i)',Pc,c);
h     = text(0.00,y,str,'FontSize',8,'FontWeight','Bold');
hPage = [hPage, h];


%-Local maxima p-values & statistics
%=======================================================================
while max(T)

	% Paginate if necessary
	%---------------------------------------------------------------
	if y < (Num+1)*dy
		h     = text(0.5,-5*dy,...
			sprintf('Page %d',spm_figure('#page')),...
			'FontSize',FS(1),'FontAngle','Italic');
		spm_figure('NewPage',[hPage,h])
		hPage = [];
		y     = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U i]   = max(T);			% largest maxima
	j       = find(A == A(i));		% maxima in cluster


    	%-Compute cluster {k} and voxel-level {u} p values for this cluster
    	%---------------------------------------------------------------
	Pz      = spm_P(1,0,   U,df,STAT,1,n);	% uncorrected p value
	Pu      = spm_P(1,0,   U,df,STAT,R,n);	% corrected     {based on Z)
	[Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n);	% [un]corrected {based on k)
	Ze      = spm_invNcdf(1 - Pz);		% Equivalent Z-variate


	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
        str   = sprintf('%-0.3f   (%-0.2f)',Pk,N(i));
	h     = text(0.18,y,str,'FontSize',8,'FontWeight','Bold');
	hPage = [hPage, h];
        str   = sprintf('%-0.3f   (%-0.2f)',Pu,Ze);
	h     = text(0.44,y,str,'FontSize',8,'FontWeight','Bold');
	hPage = [hPage, h];
        str   = sprintf('%-0.3f',Pn);
	h     = text(0.64,y,str,'FontSize',8,'FontWeight','Bold');
	hPage = [hPage, h];
        str   = sprintf('%-0.3f',Pz);
	h     = text(0.74,y,str,'FontSize',8,'FontWeight','Bold');
	hPage = [hPage, h];
        str   = sprintf('%-06.0f',M(:,i));
	h     = text(0.84,y,str,...
		'FontSize',8,'FontWeight','Bold',...
		'ButtonDownFcn','spm_mip_ui(''ShowGreens'')',...
		'Interruptible','off',...
		'UserData',M(:,i));
	hPage = [hPage, h];
 
	y     = y - dy;

	%-Print Num secondary maxima (> Dis mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-T(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
		d    = j(q(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > Dis;
		if length(D) < Num
			
			% voxel-level p values {Z}
			%-----------------------------------------------
			Pz    = spm_P(1,0,T(d),df,STAT,1,n);
			Pu    = spm_P(1,0,T(d),df,STAT,R,n);
			Ze    = spm_invNcdf(1 - Pz);

        		str   = sprintf('%-0.3f   (%-0.2f)',Pu,Ze);
			h     = text(0.44,y,str,'FontSize',8);
			hPage = [hPage, h];
        		str   = sprintf('%-0.3f',Pz);
			h     = text(0.74,y,str,'FontSize',8);
			hPage = [hPage, h];
        		str   = sprintf('%-06.0f',M(:,d));
			h     = text(0.84,y,str,'FontSize',8,...
				'ButtonDownFcn',...
				'spm_mip_ui(''ShowGreens'')',...
				'Interruptible','off',...
				'UserData',M(:,d));
			hPage = [hPage, h];
			D     = [D d];
			y     = y - dy;
		end
		end
	end
	T(j) = T(j)*0;		% Zero local maxima for this cluster
end				% end region


%-Number and register last page (if paginated)
%-----------------------------------------------------------------------
if spm_figure('#page')>1
	h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page')*[1,1]),...
		'FontSize',FS(1),'FontAngle','Italic')
	spm_figure('NewPage',[hPage,h])
end


%-Volume, resels and smoothness 
%===========================================================================
FWHMmm          = FWHM.*VOL.VOX'; 				% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n);
[P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n);

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
line([0 1],[0 0],'LineWidth',1,'Color','r')
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f (%0.3f)',u,Pz,P);
text(0.0,-1*dy,str,'FontSize',FS(1));
str = sprintf('Extent threshold {k} = %0.2f resels, p = %0.3f',k,Pn);
text(0.0,-2*dy,str,'FontSize',FS(1));
str = sprintf('Expected resels per cluster, E{n} = %0.3f',En);
text(0.0,-3*dy,str,'FontSize',FS(1));
str = sprintf('Expected number of clusters, E{m} = %0.2f',Em*Pn);
text(0.0,-4*dy,str,'FontSize',FS(1));

str = sprintf('Volume {S} = %0.0f voxels or %0.2f Resels ',VOL.S,R(length(R)));
text(0.6,-1*dy,str,'FontSize',FS(1));
str = sprintf('Degrees of freedom = %0.1f %0.1f ',df);
text(0.6,-2*dy,str,'FontSize',FS(1));
str = sprintf('Smoothness FWHM {mm} = %0.1f %0.1f %0.1f',FWHMmm);
text(0.6,-3*dy,str,'FontSize',FS(1));
str = sprintf(' {voxels} = %0.1f %0.1f %0.1f',FWHM);
text(0.6,-4*dy,str,'FontSize',FS(1));

%-End
%=======================================================================
spm('Pointer','Arrow')
