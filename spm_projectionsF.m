function  spm_projectionsF(SPMF,XYZ,u,V,Wresid,S,G,df,nInterest)
% display and analysis of SPM{F}
% FORMAT spm_projections(SPMF,XYZ,u,V,Wresid,S,G,df);
%
% SPMF- row vector of {n} voxel values
% XYZ - {3 x n} matrix of spatial locations {mm}
% u   - height threshold (uncorrected)
% V   - vector of image and voxel sizes (see spm_map.m)
% W   - vector of smoothness estimates 
%	for the original component fields 
%	(n elements for each dimension)
% S   - Lebesgue measure or volume of search space
% G   - design matrix for the experiment
% df  - degrees of freedom due to error for the uncomplete and complete
%	models.
%
%___________________________________________________________________________
%
%
% If several contrasts are specified, only voxels surviving the height
% threshold u are considered further. 
%
% The secondary maxima are selected to be at least 8mm apart.
%
%__________________________________________________________________________
% %W% Karl Friston - JBP %E%

global CWD

% If global proj_MultiPage is true then allow multi-page tables.
%---------------------------------------------------------------------------
global proj_MultiPage, if isempty(proj_MultiPage), proj_MultiPage=0; end
if proj_MultiPage, bMultiPage=1; else bMultiPage=0; end;

%-Display, characterize and print SPM{F}
%===========================================================================
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph)
spm_clf
set(Fgraph,'Units','Normalized');

%-Return if there are no voxels
%---------------------------------------------------------------------------

Q = SPMF > u;
if sum(Q) == 0
	axis off
	text(0,0.8,spm('DirTrunc',CWD));
	text(0,0.7,'No voxels above this threshold {u}','FontSize',16);
	return
end

XYZ   = XYZ(:,Q);
SPMF  = SPMF(:,Q);

%-Return if there are no voxels
%---------------------------------------------------------------------------
if sum(Q) == 0
	axis off
	text(0,0.8,spm('DirTrunc',CWD));
	text(0,0.7,'No clusters above this threshold {k}','FontSize',16);
	return
end

%-Maximium intenisty projection of SPM{F}
%===========================================================================
axes('Position',[0.05 0.5 0.5 0.5],'Tag','Empty');
spm_mip(SPMF,XYZ,V(1:6)); axis image
title('SPM{F}','FontSize',16,'Fontweight','Bold')

%-Show design matrix & parameters of interest
%===========================================================================
axes('Position',[0.65 0.6 0.2 0.2],'Tag','Empty')
imagesc(spm_DesMtxSca(G))
xlabel 'design matrix'

axes('Position',[0.65 0.8 0.2 0.1],'Tag','Empty')
C = zeros([1 size(G,2)]);
C(1:nInterest) = ones(1,nInterest);
[i j]   = bar(C);
fill(i,j,[1 1 1]*.8);
set(gca,'Xlim',[min(i) max(i)]) 
title 'Parameters of interest'; axis off; hold off

%-Display (sorted on F values and grouped by regions)
%===========================================================================

%-Characterize excursion set in terms of maxima 
%---------------------------------------------------------------------------
[N F M A] = spm_max(SPMF,XYZ,V([4 5 6]));


%-Table headings
%---------------------------------------------------------------------------
mkpt = [0 0.2 0.45 0.65];
mkpv = [0.05 0.2 0.48 0.65];

y   = 25;
axes('Position',[0.1 0.06 0.8 0.46],'Tag','Empty'); axis off
text(0,y,['P values & statistics:   ' spm('DirTrunc',CWD,40)],...
	'FontSize',12,'FontWeight','Bold');
y   = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y   = y - 1;



%-Construct tables
%---------------------------------------------------------------------------
text(mkpt(1),y,'Cluster size'       ,'FontSize',10,'Tag','Empty');
text(mkpt(2),y,'voxel-level {F}'    ,'FontSize',10,'Tag','Empty');
text(mkpt(3),y,'uncorrected'        ,'FontSize',10,'Tag','Empty');
text(mkpt(4),y,'location {mm}'      ,'FontSize',10,'Tag','Empty');
y  = y - 1;

line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y - 1;

%-Pagination variables
%-----------------------------------------------------------------------
y0 = y;
hPage = [];
nPage = 1;
Vis = 'on';


while max(F) & ((y > 2) | bMultiPage)
	% Paginate if necessary
	%-------------------------------------------------------
	if (y < 6) & proj_MultiPage
		h = text(0.5,-2,['Page ',num2str(nPage)],...
			'FontSize',8,'FontAngle','Italic',...
			'Visible',Vis);
		hPage = [hPage, h];
		spm_figure('NewPage',hPage)
		hPage = [];
		nPage = nPage + 1;
		Vis   = 'off';
		y = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U i] = max(F);				% largest F value
	j     = find(A == A(i));		% maxima in cluster

    	%-Compute cluster {k} and voxel-level p values for this cluster
    	%---------------------------------------------------------------
	Pu  = spm_pF(S, Wresid, df, U);		% voxel-level p value
	Pz  = 1 - spm_Fcdf(U,df);		% uncorrected p value

%	Pkn   = spm_Pkn(N(i),U-u,W,u,S);	% Bivariate cluster-height p
%	Pk    = spm_P(1,W,u,N(i),S);		% cluster-level p value
%	Pu    = spm_P(1,W,U,0,S);		% voxel-level p value
%	Pz    = 1 - spm_Ncdf(U);		% uncorrected p value

	%-Print cluster and maximum voxel-level p values {F}
    	%---------------------------------------------------------------

        str   = sprintf('%-2.0i',N(i));
	h = text(mkpv(1),y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f   (%-0.2f)',Pu,U);
	h = text(mkpv(2),y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-0.3f',Pz);
	h = text(mkpv(3),y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
        str   = sprintf('%-6.0f',M(:,i));
	h = text(mkpv(4),y,str,'FontSize',8,'FontWeight','Bold','Visible',Vis);
	hPage = [hPage, h];
 
	y     = y - 1;

	%-Print 2 secondary maxima (> 8mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-F(j));				% sort on F value
	D     = i;
	for i = 1:length(q)
		d    = j(q(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > 8;
		    if length(D) < 3
			
			% voxel-level p values {F}
			%-----------------------------------------------

			Pu  = spm_pF(S, Wresid, df, F(d));% voxel-level p value
			Pz  = 1 - spm_Fcdf(F(d),df);	% uncorrected p value

        		str = sprintf('%-0.3f   (%-0.2f)',Pu,F(d));
			h = text(mkpv(2),y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
        		str = sprintf('%-0.3f',Pz);
			h = text(mkpv(3),y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
        		str = sprintf('%-6.0f',M(:,d));
			h = text(mkpv(4),y,str,'FontSize',8,'Visible',Vis);
			hPage = [hPage, h];
			D   = [D d];
			y   = y - 1;
		    end
		end
	end
	F(j) = F(j)*0;			% Zero local maxima for this cluster
end					% end region

%-Number and paginate last page (if pagination was required)
if strcmp(Vis,'off')
	%-Label last page
	h = text(0.5,-2,['Page ',num2str(nPage)],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis);
	hPage = [hPage, h];
	spm_figure('NewPage',hPage);
end

y = min(3,y);

%-Footnote with SPM parameters
%-----------------------------------------------------------------------

%PWD		= pwd;
%PWD    	= PWD([max(find(pwd == '/')):length(pwd)]);
D           	= length(mean(Wresid));				% dimension
FWHM        	= sqrt(8*log(2))*mean(Wresid).*V(([1:D] + 3),1)'; % FWHM {mm}
RESEL      	= S*prod(V([1:D] + 3))/prod(FWHM);		% RESELS



y  = y - .5;
line([0 1],[y y],'LineWidth',1,'Color',[0 0 0]);

y  = y - .5;
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f',u,1 - spm_Fcdf(u,df));
text(0,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Degrees of freedom  = %0.0i ; %0.0i ',df(1),df(2));
text(0,y,str,'FontSize',8);
y  = y - 1;
if D == 2
    str = sprintf('Smoothness = %0.1f %0.1f mm {FWHM}',FWHM);
end
if D == 3
    str = sprintf('Smoothness = %0.1f %0.1f %0.1f mm {FWHM}',FWHM);
end
text(0,y,str,'FontSize',8);

y  = y - .5;
line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])

