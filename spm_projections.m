function spm_projections(Z,df,STAT,R,S,FWHM,u,k,XYZ,V,G,C)
% display and analysis of SPM{Z} [SPM{t}]
% FORMAT spm_projections(Z,df,STAT,R,S,FWHM,u,k,XYZ,V,G,C)
%
% Z    - matrix of {n x m} statistic voxel values
% df   - [df{interest} df{error}]
% STAT - statistic ('Z', 'T', 'X' or 'F')
% R    - RESEL count defining search volume {resels}
% S    - Lebesgue measure or volume of search space {voxels}
% FWHM - vector of effective smoothness estimates {voxels}
% u    - height threshold
% k    - extent threshold {resels}
% XYZ  - {3 x m} matrix of spatial locations {mm}
% V    - vector of image and voxel sizes (see spm_map.m)
% G    - design matrix for the experiment
% C    - row vector contrast[s]
%
%___________________________________________________________________________
%
% spm_projections characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected subsets
% or clusters).  The significance of the results are based on set, cluster
% and voxel-level inferences using distributional approximations from the
% Theory of Gaussian Feilds.  These distributions assume that the SPM is
% a reasonable lattice approximation with known component field smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, voxels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C = P(u,k,C).  For each cluster
% of size K the cluster-level inference is based on P(u,K,1) and for each
% voxel (or selected maxima) of height U, in that cluster, the voxel-level
% inference is based on P(U,0,1).  All three levels of inference are
% supported with a tabular presentation of the p values and the underlying
% statistic (u, k or c).  The table is grouped by regions and sorted on
% the Z-variate of the primary maxima.  See 'Sections' in the help facility
% and spm_maxima.m for a more complete characterization of maxima within
% a region. Z-variates (based on the uncorredcted p value) are the Z score
% equivalent of the statistic and are preferred for tabular reporting.
%
% The secondary maxima are selected to be at least 8mm apart.
%
%______________________________________________________________________
% %W% Karl Friston %E%

global CWD

%-If global proj_MultiPage is true then allow multi-page tables.
%-----------------------------------------------------------------------
global proj_MultiPage, if isempty(proj_MultiPage), proj_MultiPage = 0; end


%-Display, characterize and print SPM{Z}
%=======================================================================

% Get figure handles
%-----------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph)
spm_clf
set(Fgraph,'Units','Normalized');

%-find minima of SPM
%-----------------------------------------------------------------------
n      = size(Z,1);
Z      = min(Z,[],1);

%-Maximium intensity projection of SPM
%-----------------------------------------------------------------------
axes('Position',[0.05 0.5 0.5 0.5],'Tag','Empty');
spm_mip(Z,XYZ,V); axis image
title(['SPM{' STAT '}'],'FontSize',16,'Fontweight','Bold')


%-Show design matrix & contrast
%-----------------------------------------------------------------------
axes('Position',[0.65 0.6 0.2 0.2],'Tag','Empty')
imagesc(spm_DesMtx('sca',G))
xlabel 'design matrix'
dy    = 0.1/size(C,1);
for i = 1:size(C,1)
	axes('Position',[0.65 (0.8 + dy*(i - 1)) 0.2 dy],'Tag','Empty')
	bar(C(i,:));
	axis tight, axis off
end



%-Display (sorted on Z values and grouped by regions)
%=======================================================================

%-Characterize excursion set in terms of maxima 
%-----------------------------------------------------------------------
[N T M A] = spm_max(Z,XYZ,V([4 5 6]));
N         = N*max(R)/S; 			      % voxels -> resels


%-Table headings
%-----------------------------------------------------------------------
y     = 25;
axes('Position',[0.1 0.06 0.8 0.46],'Tag','Empty'); axis off
text(0,y,['P values & statistics:   ' spm_str_manip(CWD,'a40')],...
	'FontSize',12,'FontWeight','Bold');
y     = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y     = y - 1;

%-Construct tables
%-----------------------------------------------------------------------
text(0.00,y,'set-level {c}'      ,'FontSize',10,'Tag','Empty');
text(0.18,y,'cluster-level {k}  ','FontSize',10,'Tag','Empty');
text(0.42,y,'voxel-level {Z}'    ,'FontSize',10,'Tag','Empty');
text(0.62,y,'uncorrected k & Z'  ,'FontSize',10,'Tag','Empty');
text(0.86,y,'x,y,z {mm}'         ,'FontSize',10,'Tag','Empty');


%-Pagination variables
%-----------------------------------------------------------------------
y     = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y     = y - 1; y0 = y; hPage = []; nPage = 1; Vis = 'on';

%-Set-level p values {c}
%-----------------------------------------------------------------------
c     = max(A);					% number of clusters
Pc    = spm_P(c,k,u,df,STAT,R,n);		% set-level p value

str   = sprintf('%-0.3f   (%i)',Pc,c);
h     = text(0.00,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
hPage = [hPage, h];

while max(T) & ((y > 2) | proj_MultiPage)

	% Paginate if necessary
	%---------------------------------------------------------------
	if (y < 6) & proj_MultiPage
		h     = text(0.5,-2,['Page ',num2str(nPage)],...
			'FontSize',8,'FontAngle','Italic',...
			'Visible',Vis);
		hPage = [hPage, h];
		spm_figure('NewPage',hPage)
		hPage = [];
		nPage = nPage + 1;
		Vis   = 'off';
		y     = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U i]   = max(T);			% largest maxima
	j       = find(A == A(i));		% maxima in cluster

    	%-Compute cluster {k} and voxel-level {u} p values for this cluster
    	%---------------------------------------------------------------
	Pz      = spm_P(1,0,   U,df,STAT,1,n);	% uncorrected p value
	Pu      = spm_P(1,0,   U,df,STAT,R,n);	% corrected   {based on height)
	[Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n);	% [un]corrected {based on k)
	Ze      = 1 - spm_invNcdf(Pz);		% Equivalent Z-variate

        if N(i) < k; break; end

	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
        str   = sprintf('%-0.3f   (%-0.2f)',Pk,N(i));
	h     = text(0.18,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f   (%-0.2f)',Pu,Ze);
	h     = text(0.44,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f',Pn);
	h     = text(0.64,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f',Pz);
	h     = text(0.74,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-6.0f',M(:,i));
	h     = text(0.84,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
 
	y     = y - 1;

	%-Print 2 secondary maxima (> 8mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-T(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
		d    = j(q(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > 8;
		if length(D) < 3
			
			% voxel-level p values {Z}
			%-----------------------------------------------
			Pz    = spm_P(1,0,T(d),df,STAT,1,n);
			Pu    = spm_P(1,0,T(d),df,STAT,R,n);
			Ze    = 1 - spm_invNcdf(Pz);

        		str   = sprintf('%-0.3f   (%-0.2f)',Pu,Ze);
			h     = text(0.44,y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
        		str   = sprintf('%-0.3f',Pz);
			h     = text(0.74,y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
        		str   = sprintf('%-6.0f',M(:,d));
			h     = text(0.84,y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
			D     = [D d];
			y     = y - 1;
		end
		end
	end
	T(j) = T(j)*0;		% Zero local maxima for this cluster
end					% end region

%-Number and paginate last page (if pagination was required)
%---------------------------------------------------------------------------
if strcmp(Vis,'off')
	%-Label last page
	h = text(0.5,-2,['Page ',num2str(nPage)],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis);
	hPage = [hPage, h];
	spm_figure('NewPage',hPage);
end

y      = min(3,y);

%-Volume, resels and smoothness 
%===========================================================================
FWHMmm          = FWHM.*V(([1:length(FWHM)] + 3),1)'; 		% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n);
[P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n)

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
y   = y + 0.5;
line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])
y   = y - 0.5;
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f (%0.3f)',u,Pz,P);
text(0,y,str,'FontSize',8);
y   = y - 1;
str = sprintf('Extent threshold {k} = %0.2f resels, p = %0.3f',k,Pn);
text(0,y,str,'FontSize',8);
y   = y - 1;
str = sprintf('Expected resels per cluster, E{n} = %0.1f',En);
text(0,y,str,'FontSize',8);
y   = y - 1;
str = sprintf('Expected number of clusters, E{m} = %0.1f',Em*Pn);
text(0,y,str,'FontSize',8);
y   = y + 3;
str = sprintf('Volume {S} = %i voxels or %0.1f Resels ',S,max(R));
text(0.6,y,str,'FontSize',8);
y   = y - 1;
str = sprintf('Degrees of freedom = %0.1f %0.1f ',df);
text(0.6,y,str,'FontSize',8);
y   = y - 1;
str = sprintf('Smoothness FWHM {mm} = %0.1f %0.1f %0.1f',FWHMmm);
text(0.6,y,str,'FontSize',8);
y   = y - 1;
str = sprintf(' {voxels} = %0.1f %0.1f %0.1f',FWHM);
text(0.6,y,str,'FontSize',8);

