function [Z,XYZ] = spm_projections(Z,XYZ,u,k,V,W,S,G,C,df)
% display and analysis of SPM{Z} [SPM{t}]
% FORMAT [Z,XYZ] = spm_projections(Z,XYZ,u,k,V,W,S,G,C,df);
%
% Z   - row vector of {n} Z or F statistic voxel values
% XYZ - {3 x n} matrix of spatial locations {mm}
% u   - height threshold (uncorrected)
% k   - extent threshold (uncorrected)
% V   - vector of image and voxel sizes (see spm_map.m)
% W   - vector of smoothness estimates (one element for each dimension)
% S   - Lebesgue measure or volume of search space
% G   - design matrix for the experiment
% C   - row vector contrast for the SPM{Z}
% df  - degrees of freedom due to error 
% NB:   if df = [df1 df2] then F statistics are assumed (otherwise Z).
%
% Z & XYZ - filtered according to threshold criteria
%___________________________________________________________________________
%
% spm_projections applies thresholds {u & k} to a point list of voxel values
% (specified with their locations {XYZ}) and characterizes the resulting
% excursion set as a collection of face, edge and vertex connected subsets
% or clusters.  The significance of the results are based on set, cluster
% and voxel-level inferences using distributional approximations from the
% Theory of Gaussian Feilds.  These distributions assume that the SPM is
% a reasonable lattice approximation to a Gaussian field of smoothness
% {W = FWHM/sqrt(8.ln2)} (or in the instance of SPM{F} the component feilds
% have this smoothness
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
% the Z value of the primary maxima.  See 'Sections' in the help facility
% and spm_maxima.m for a more complete characterization of maxima within
% a region. For SPM{F} only voxel-level inferences are available at
% the present time. NaN means not a number.  When expressions for these
% corrected p values become available the NaNs will be replaced by
% proper values.
%
% The secondary maxima are selected to be at least 8mm apart.
%
%______________________________________________________________________
% %W% Karl Friston %E%

global CWD

%-If called with output arguments (i.e. from spm_projections_ui.m)
% *and* global proj_MultiPage is true then allow multi-page tables.
%-----------------------------------------------------------------------
global proj_MultiPage, if isempty(proj_MultiPage), proj_MultiPage = 0; end
if (nargout & proj_MultiPage), bMultiPage = 1; else, bMultiPage = 0;   end


%-Display, characterize and print SPM{Z}
%=======================================================================

% Get figure handles
%-----------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph)
spm_clf
set(Fgraph,'Units','Normalized');

% SPM{Z} or SPM{F}
%-----------------------------------------------------------------------
SPMZ = length(df) == 1;
SPMF = length(df) == 2;

% eliminate voxels on the basis of u
%-----------------------------------------------------------------------
Q     = find(Z > u);
if ~length(Q); return; end
Z     = Z(Q);
XYZ   = XYZ(:,Q);

% eliminate voxels on the basis of k
%-----------------------------------------------------------------------
A     = spm_clusters(XYZ,V([4 5 6]));
Q     = [];
for i = 1:max(A)
	j = find(A == i);
	if length(j) >= k
		Q = [Q j];
	end
end
if ~length(Q); return; end
Z     = Z(Q);
XYZ   = XYZ(:,Q);


%-Maximium intensity projection of SPM{Z}
%-----------------------------------------------------------------------
axes('Position',[0.05 0.5 0.5 0.5],'Tag','Empty');
spm_mip(Z,XYZ,V(1:6)); axis image
if SPMZ
	title('SPM{Z}','FontSize',16,'Fontweight','Bold')
elseif SPMF
	title('SPM{F}','FontSize',16,'Fontweight','Bold')
end


%-Show design matrix & contrast
%-----------------------------------------------------------------------
axes('Position',[0.65 0.6 0.2 0.2],'Tag','Empty')
imagesc(spm_DesMtxSca(G))
xlabel 'design matrix'
if SPMZ
	dy    = 0.1/size(C,1);
	for i = 1:size(C,1)
		axes('Position',[0.65 (0.8 + dy*(i - 1)) 0.2 dy],'Tag','Empty')
		[p q]   = bar(C(i,:));
		fill(p,q,[1 1 1]*.8);
		set(gca,'Xlim',[min(p) max(p)]); axis off
	end
	title 'contrast'
end


%-Display (sorted on Z values and grouped by regions)
%=======================================================================

%-Characterize excursion set in terms of maxima 
%-----------------------------------------------------------------------
[N Zm M A] = spm_max(Z,XYZ,V([4 5 6]));


%-Table headings
%-----------------------------------------------------------------------
y   = 25;
axes('Position',[0.1 0.06 0.8 0.46],'Tag','Empty'); axis off
text(0,y,['P values & statistics:   ' spm('DirTrunc',CWD,40)],...
	'FontSize',12,'FontWeight','Bold');
y   = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y   = y - 1;

%-Construct tables
%-----------------------------------------------------------------------
if SPMF
	text(0.00,y,'set-level {c}'      ,'FontSize',10,'Tag','Empty');
	text(0.18,y,'cluster-level {k,F}','FontSize',10,'Tag','Empty');
	text(0.42,y,'voxel-level {F}'    ,'FontSize',10,'Tag','Empty');
	text(0.62,y,'uncorrected k & F'  ,'FontSize',10,'Tag','Empty');
	text(0.86,y,'x,y,z {mm}'         ,'FontSize',10,'Tag','Empty');

elseif SPMZ

	text(0.00,y,'set-level {c}'      ,'FontSize',10,'Tag','Empty');
	text(0.18,y,'cluster-level {k,Z}','FontSize',10,'Tag','Empty');
	text(0.42,y,'voxel-level {Z}'    ,'FontSize',10,'Tag','Empty');
	text(0.62,y,'uncorrected k & Z'  ,'FontSize',10,'Tag','Empty');
	text(0.86,y,'x,y,z {mm}'         ,'FontSize',10,'Tag','Empty');
end

%-Pagination variables
%-----------------------------------------------------------------------
y     = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y     = y - 1;
y0    = y;
hPage = [];
nPage = 1;
Vis   = 'on';

%-Set-level p values {c}
%-----------------------------------------------------------------------
c     = max(A);					% number of clusters
if SPMF
	Pc = NaN;				% set-level p value

elseif SPMZ

	Pc = spm_P(c,W,u,k,S);			% set-level p value

end
str   = sprintf('%-0.3f   (%i)',Pc,c);
h     = text(0.00,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
hPage = [hPage, h];

while max(Zm) & ((y > 2) | bMultiPage)

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
	[U i] = max(Zm);			% largest Z maxima
	j     = find(A == A(i));		% maxima in cluster

        if SPMZ

    	%-Compute cluster {k} and voxel-level p values for this cluster
    	%---------------------------------------------------------------
	Pk    = spm_P(1,W,u,N(i),S);		% cluster-level p value
	Pkn   = spm_Pkn(N(i),(U - u),W,u,S);	% cluster-level (Bivariate)
	Pu    = spm_P(1,W,U,0,S);		% voxel-level p value
	Pn    = 1 - spm_kcdf(N(i),u,W);		% uncorrected p value (k)
	Pz    = 1 - spm_Ncdf(U);		% uncorrected p value (Z)

	elseif SPMF

    	%-Compute cluster {k} and voxel-level p values for this cluster
    	%---------------------------------------------------------------
	Pk    = NaN;				% cluster-level p value
	Pkn   = NaN;				% cluster-level (Bivariate)
	Pu    = spm_pF(S,W,df,U);		% voxel-level p value
	Pn    = NaN;				% uncorrected p value (k)
	Pz    = 1 - spm_Fcdf(U,df);		% uncorrected p value (F)

	end

        if N(i) < k; break; end

	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
        str   = sprintf('%-0.3f   (%i, %0.2f)',Pkn,N(i),U);
	h     = text(0.18,y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f   (%-0.2f)',Pu,U);
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
	[l q] = sort(-Zm(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
		d    = j(q(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > 8;
		    if length(D) < 3
			
			% voxel-level p values {Z}
			%-----------------------------------------------
			if SPMZ

			Pu    = spm_P(1,W,Zm(d),0,S);	% voxel-level p value
			Pz    = 1 - spm_Ncdf(Zm(d));	% uncorrected p value

			elseif SPMF

			Pu    = spm_pF(S,W,df,Zm(d));	% voxel-level p value
			Pz    = 1 - spm_Fcdf(Zm(d),df);	% uncorrected p value

			end

        		str   = sprintf('%-0.3f   (%-0.2f)',Pu,Zm(d));
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
	Zm(j) = Zm(j)*0;		% Zero local maxima for this cluster
end					% end region

%-Number and paginate last page (if pagination was required)
%-----------------------------------------------------------------------
if strcmp(Vis,'off')
	%-Label last page
	h = text(0.5,-2,['Page ',num2str(nPage)],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis);
	hPage = [hPage, h];
	spm_figure('NewPage',hPage);
end

y = min(3,y);

%-Volume, resels and smoothness 
%=======================================================================
D               = length(W);					% dimension
FWHM            = sqrt(8*log(2))*W.*V(([1:D] + 3),1)'; 		% FWHM {mm}
RESEL           = S*prod(V([1:D] + 3))/prod(FWHM);		% RESELS

if SPMZ

	Pu              = 1 - spm_Ncdf(u);
	Pn              = 1 - spm_kcdf(k,u,W);
	[P,EN,Em,En,Pk] = spm_P(1,W,u,k,S);

elseif SPMF

	Pu = 1 - spm_Fcdf(u,df);
	Pn = NaN;
	Em = -log( 1 - spm_pF(S,W,df,u));
	En = NaN;
	Pk = 1;
end


%-Footnote with SPM parameters
%-----------------------------------------------------------------------
y  = y + 0.5;
line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])
y  = y - 0.5;
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f',u,Pu);
text(0,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Extent threshold {k} = %i voxels, p = %0.3f',k,Pn);
text(0,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Expected voxels per cluster, E{n} = %0.1f',En);
text(0,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Expected number of clusters, E{m} = %0.1f',Em*Pk);
text(0,y,str,'FontSize',8);
y  = y + 3;
str = sprintf('Volume {S} = %i voxels or %0.1f Resels ',S,RESEL);
text(0.6,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Degrees of freedom due to error = %i %i ',df);
text(0.6,y,str,'FontSize',8);
y  = y - 1;

if D == 2
    str = sprintf('Smoothness = %0.1f %0.1f mm {FWHM}',FWHM);
end
if D == 3
    str = sprintf('Smoothness = %0.1f %0.1f %0.1f mm {FWHM}',FWHM);
end
text(0.6,y,str,'FontSize',8);
y  = y - 1;
if D == 2
    str = sprintf(' = %0.1f %0.1f {voxels}',W.*V([1:D] + 3)');
end
if D == 3
    str = sprintf(' = %0.1f %0.1f %0.1f {voxels}',W.*V([1:D] + 3)');
end
text(0.6,y,str,'FontSize',8);

