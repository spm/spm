function [t,XYZ]=spm_projections(t,XYZ,U,V,W,S,G,C,df,pV,nam)
% display and analysis of SPM{Z} [SPM{t}]
% FORMAT spm_projections(Z,XYZ,U,V,W,S,G,C,df,[pV,nam])
% t   - row vector of {n} voxel values
% XYZ - {3 x n} matrix of spatial locations {mm}
% U   - threshold (uncorrected)
% V   - vector of image and voxel sizes (see spm_map.m)
% W   - vector of smoothness estimates (one element for each dimension)
% S   - Lebesgue measure or volume of search space
% G   - design Matrix for the experiment
% C   - row vector contrast for the SPM{Z}
% df  - degrees of freedom due to error
% pV  - threshold (corrected)
% nam - design matrix labels
%-t & XYZ are filtered by corrected p-value, and returned
%___________________________________________________________________________
%
% spm_projections applies a threshold {U} to a point list of voxel values {Z}
% (specified with their locations {XYZ}) and characterizes the resulting
% excursion set as a collection of face, edge and vertex connected subsets
% or regions.  The significance of each region is assessed on the basis of its
% spatial extent {k} [P(nmax > k)] and its peak Z value {u} [P(Zmax > u)].  
% These estimates are based on distributional approximations assuming that the
% voxel values come from a Gaussian field of smoothness {W = FWHM/sqrt(8.ln2)}.
%
% The approximation for P(Zmax > u) is based upon <m>, the expected
% number of maxima m  above u.  Assuming a Poisson form for the p.d.f of the
% number of maxima the probability of there being one or more maxima (i.e.
% 1 - P(m = 0) = 1 - exp(-<m>) is calculated using expressions for <m> based
% on the theory of Gaussian Fields.  An extended analysis, using 
% asymptotically correct distributional forms and <m> gives the equivalent
% probability of one or more regions of the size observed, or larger.
% These probabilities can be regarded as corrected p values.
%
% For each region the size {k} and significance [P(nmax > k)] are tabulated
% The primary maximum (largest Z value) and up to 2 secondary maxima are
% provided with their corrected [P(Zmax > u)] and uncorrected p values.
% The table is grouped by regions and sorted on the Z value of the primary 
% maxima.  See 'Sections' in the help facility and spm_maxima.m for a more
% complete characterization of maxima within a region.
%
% If several contrasts are specified the uncorrected p value is raised
% to the appropriate power, assuming the contrasts are othogonal.  Corrected
% p values are not, in this instance, shown.
%
% The secondary maxima are selected to be at least 8mm apart.
%
% The uncorrected threshold defines the excursion set and the corrected
% threshold removes regions that fail to reach a corrected level of
% significance.  To see all the voxels set the corrected p value to 1.
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
global CWD

% display, characterize and print SPM{t}
%===========================================================================
figure(3); spm_clf; set(3,'Units','Normalized');

% apply threshold {U}
%---------------------------------------------------------------------------
if size(t,1) > 1
	Q = all(t > U);
else
	Q = t > U;
end
t       = t(1,Q);
XYZ     = XYZ(:,Q);

% return if there are no voxels
%---------------------------------------------------------------------------
if sum(Q) == 0
	axis off
	text(0,0.8,CWD,'Fontsize',16,'FontWeight','Bold');
	text(0,0.7,'No voxels above this threshold');
	return
end

% set the corrected p value and design matrix names (if not specified)
%---------------------------------------------------------------------------
if nargin < 10; pV  = 1;  end
if nargin < 11; nam = []; end


% filter on P(nmax > k) and P(Zmax > u)
%---------------------------------------------------------------------------
A         = spm_clusters(XYZ,V([4 5 6]));
Q         = [];
for i     = 1:max(A)
	j = find(A == i);
	if ( spm_Pn(length(j),W,U,S) <= pV  ) | (spm_Pz(W,max(t(j)),S) <= pV )
		Q = [Q j];
	end
end

% return if there are no voxels
%---------------------------------------------------------------------------
if sum(Q) == 0
	axis off
	text(0,0.8,CWD,'Fontsize',16,'FontWeight','Bold');
	text(0,0.7,'No voxels significant at this level');
	return
end

t       = t(1,Q);
XYZ     = XYZ(:,Q);



% maximium intenisty projection of SPM{Z}
%---------------------------------------------------------------------------
axes('Position',[0.05 0.5 0.5 0.5])
spm_mip(t,XYZ,V(1:6)); axis image
title('SPM{Z}','FontSize',16,'Fontweight','Bold')

axes('Position',[0.65 0.6 0.2 0.2])
imagesc((spm_DesMtxSca(G,nam) + 1)*32)
xlabel 'Design Matrix'
axes('Position',[0.65 0.8 0.2 0.1])
bar(C'); hold on
[j k]   = bar(C(1,:));
fill(j,k,[1 1 1]*.8);
set(gca,'Xlim',[min(j) max(j)]) 
title 'contrast'; axis off; hold off

	
% characterize local excursions in terms of maxima {N - voxels},  maxima
% Z values {Z} and locations (M)
%---------------------------------------------------------------------------
[N Z M A] = spm_max(t,XYZ,V([4 5 6]));
Pn        = spm_Pn(N,W,U,S);
Pz        = spm_Pz(W,Z,S);
Pu        = 1 - spm_Ncdf(Z);

% display (sorted on Z values and grouped by regions)
%===========================================================================


% table headings
%---------------------------------------------------------------------------
axes('Position',[0.1 0.06 0.8 0.46]); axis off
y  = 24;
text(0,y,['Regional effects:    ' CWD],'Fontsize',16,'FontWeight','Bold');
y  = y - 1.2;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y - 1;



% multiple orthogonal contrasts
%===========================================================================
if (size(C,1) > 1) & all(all(~(C*C' - diag(diag(C*C')))))

    % product of p values
    %-----------------------------------------------------------------------
    Pu  = Pu*( (1 - spm_Ncdf(U))^(size(C,1) - 1) );

    text(0.00,y,'region');
    text(0.10,y,'size {k}');
    text(0.42,y,'Z');
    text(0.54,y,'P value (Uncorrected)','Fontsize',10);
    text(0.84,y,'{x,y,z mm}');
    y  = y - 1;

    line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
    y  = y - 1;

    % list of maxima
    %-----------------------------------------------------------------------
    region = 1;
    while max(Z) & (y > 0)

	[j i] = max(Z);				% largest Z value
	j     = find(A == A(i));		% maxima in that region


	% print region and largest maximum
	%-------------------------------------------------------------------
	text(0.00,y,sprintf('%-0.0f',region),'Fontsize',10,'FontWeight','Bold')
	text(0.10,y,sprintf('%-0.0f',N(i))  ,'Fontsize',10,'FontWeight','Bold')
	text(0.42,y,sprintf('%-0.2f',Z(i))  ,'Fontsize',10,'FontWeight','Bold')
	text(0.82,y,sprintf('%-6.0f',M(:,i)),'Fontsize',10,'FontWeight','Bold')
	text(0.60,y,sprintf('%-0.3g',Pu(i)) ,'Fontsize',10)
	y     = y - 1;

	% print region and 4 secondary maxima (8mm apart)
	%-------------------------------------------------------------------
	[l k] = sort(-Z(j));			% sort on Z value
	D     = i;
	for i = 1:length(k)
		d    = j(k(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > 8;
		    if length(D) < 3
			text(0.42,y,sprintf('%-0.2f',Z(d))   ,'Fontsize',10)
			text(0.82,y,sprintf('%-6.0f',M(:,d)) ,'Fontsize',10)
			text(0.60,y,sprintf('%-0.3g',Pu(d))  ,'Fontsize',10)
			D = [D d];
			y = y - 1;
		    end
		end
	end
	Z(j) = Z(j)*0;
	region = region + 1;			% next region
    end						% end region
    text(0,y,'NB multiple orthogonal contrasts','Color',[1 0 0])
    y = y - 1;


else			% single contrast or nonorthogonal contrasts
%===========================================================================
    text(0.00,y,'region');
    text(0.10,y,'size {k}');
    text(0.24,y,'P(n      > k)');
    text(0.29,y,'max','Fontsize',8);
    text(0.42,y,'Z');
    text(0.48,y,'P(Z      > u)');
    text(0.53,y,'max','Fontsize',8);
    text(0.64,y,'(Uncorrected)','Fontsize',10);
    text(0.84,y,'{x,y,z mm}');
    y  = y - 1;

    line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
    y  = y - 1;

    % list of maxima
    %-----------------------------------------------------------------------
    region = 1;
    while max(Z) & (y > 0)

	[j i] = max(Z);				% largest Z value
	j     = find(A == A(i));		% maxima in that region


	% print region and largest maximum
	%-------------------------------------------------------------------
	text(0.00,y,sprintf('%-0.0f',region),'Fontsize',10,'FontWeight','Bold')
	text(0.10,y,sprintf('%-0.0f',N(i))  ,'Fontsize',10,'FontWeight','Bold')
	text(0.24,y,sprintf('%-0.3f',Pn(i)) ,'Fontsize',10,'FontWeight','Bold')
	text(0.42,y,sprintf('%-0.2f',Z(i))  ,'Fontsize',10,'FontWeight','Bold')
	text(0.54,y,sprintf('%-0.3f',Pz(i)) ,'Fontsize',10,'FontWeight','Bold')
	text(0.82,y,sprintf('%-6.0f',M(:,i)),'Fontsize',10,'FontWeight','Bold')
	text(0.64,y,sprintf('(%-0.3f)',Pu(i)) ,'Fontsize',10)
	y     = y - 1;

	% print region and 4 secondary maxima (8mm apart)
	%-------------------------------------------------------------------
	[l k] = sort(-Z(j));			% sort on Z value
	D     = i;
	for i = 1:length(k)
		d    = j(k(i));
		if min(sqrt(sum((M(:,D) - M(:,d)*ones(1,size(D,2))).^2))) > 8;
		    if length(D) < 3
			text(0.42,y,sprintf('%-0.2f',Z(d))   ,'Fontsize',10)
			text(0.54,y,sprintf('%-0.3f',Pz(d))  ,'Fontsize',10)
			text(0.82,y,sprintf('%-6.0f',M(:,d)) ,'Fontsize',10)
			text(0.64,y,sprintf('(%-0.3f)',Pu(d)),'Fontsize',10)
			D = [D d];
			y = y - 1;
		    end
		end
	end
	Z(j) = Z(j)*0;
	region = region + 1;			% next region
    end						% end region
end						% end conditional


line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])
y     = y - 1;

% volume, resels and smoothness 
%---------------------------------------------------------------------------
D     = length(W);					% dimension of SPM
FWHM  = sqrt(8*log(2))*W.*V(([1:D] + 3),1)'; 		% FWHM in mm
RESEL = S*prod(V([1:D] + 3))/prod(FWHM);		% RESELS

% footnote with SPM parameters
%---------------------------------------------------------------------------
d  = sprintf('Threshold = %0.2f; Volume [S] = %0.0f voxels; df = %0.0f',U,S,df);
text(0,y,d,'Fontsize',10);
y  = y - 1;
if D == 2
    d = sprintf('FWHM = [%0.1f %0.1f] mm (i.e. %0.0f RESELS) ',FWHM,RESEL);
end
if D == 3
    d = sprintf('FWHM = [%0.1f %0.1f %0.1f] mm (i.e. %0.0f RESELS)',FWHM,RESEL);
end
text(0,y,d,'Fontsize',10);


set(gca,'Ylim',[y 24]);
