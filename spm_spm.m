function spm_spm(V,K,H,C,B,G,CONTRAST,Ut,ORIGIN,TH,FLIP,names);
% Statistical analysis with the General linear model
% FORMAT spm_spm(V,K,H,C,B,G,CONTRAST,Ut,ORIGIN,TH,FLIP,names);
% V   - {12 x q} matrix of identifiers of memory mapped data {q scans}
%
% K   - {q  x 1} constant  subpartition of the design matrix {if it exists}
% H   - {q  x h} condition subpartition of the design matrix {h conditions}
% C   - {q  x c} covariate subpartition of the design matrix {c covariates}
% B   - {q  x n} block     subpartition of the design matrix {n subjects}
% G   - {q  x g} confound  subpartition of the design matrix {g covariates}
%
% names      - string matrix holding the names of the effects
% CONTRAST   - matrix of contrasts, one per row, with p elements.
% Ut         - threshold for SPM{Z}
% ORIGIN     - the voxel correpsonding to [0 0 0] in mm
% TH         - thresholds for each image defining voxels of interest
% FLIP       - Left-right orientation
%	           - 0 left = left  - neurological convention
%	           - 1 left = right - radiological convention
%_______________________________________________________________________
%
% spm_spm is the heart of the SPM package and implements the general
% linear model in terms of a design matrix (composed of K H C B
% and G) and the data (V).  Significant compounds of the estimated
% parameters are assessed with a quotient that has the t distribution
% under the null hyypothesis.  The resulting SPM{t} is transformed to
% the Unit Gaussian distribution [SPM{Z}] and characterized by further
% analysis using the theory of Gaussian Fields (see spm_projections.m
% for more details)
%
% The outputs of this routine are a series of .mat files containing
% paramter estimates, adjusted values, SPM{Z} etc that are written to
% CWD (see spm_defults.m).  IMPORTANT: Existing results are overwritten
% without prompting
%
% Voxels are retained for further analysis if the F ratio for that
% voxel is significant (p < UFp uncorrected) and all the voxels have a
% reasonably high activity [the threhsold is specified as a fraction
% (usually 0.8) of the whole brain mean].
%
%   SPMF.mat contains a 1 x N vector of F values reflecting the omnibus
% significance of effects [of interest] at each of the N 'significant'
% voxels.  'Significance' is defined by the p-value of the F threshold
% (p < UFp).
%
%   XYZ.mat contains a 3 x N matrix of the x,y and z location of the
% voxels in SPMF in mm (usually referring the the standard anatomical
% space (Talairach and Tournoux 1988)} (0,0,0) corresponds to the
% centre of the voxel specified by ORIGIN in the *.hdr of the original
% and related data.
%
%   BETA.mat contains a p x N matrix of the p parameter estimates at
% each of the N voxels.  These parameters include all effects
% specified by the design matrix.
%
%   XA.mat contains a q x N matrix of adjusted activity values having 
% removed the effects of no interest at each of the N voxels for all q
% scans.
%
%   SPMt.mat contains a c x N matrix of the c SPM{Z} defined by the c  
% contrasts supplied for all N voxels at locations XYZ.
%
%   SPM.mat  contains a collection of matrices that pertain to the 
% analysis; including the partitions of the design matrix [K H C B G], the
% number of voxels analyzed (S), image and voxel dimensions [V],
% smoothness estimates of the SPM{Z} [W FWHM], threshold for SPM{F}
% [UF] and the contrasts used [CONTRAST].  See below for a complete
% listing.
%
% Output to the results window includes maximum intensity projections
% of the SPM{F}, the design matrix and a series of pages for the SPM{Z}
% (see 'Results' in the help application).
%
% Variables saved in SPM.mat
%-----------------------------------------------------------------------
% K	-	column of ones {constant partition}
% H	-	condition partition of design matrix
% C	-	covariate partition of design matrix
% B	-	block     partition of design matrix
% G	-	confound  partition of design matrix
% S	-	Lebegue measure or volume {voxels}
% UF	-	Threshold for F ratio of variances
% V(1)	-	x image size {voxels}
% V(2)	-	y image size {voxels}
% V(3)	-	z image size {voxels}
% V(4)	-	x voxel size {mm}
% V(5)	-	y voxel size {mm}
% V(6)	-	z voxel size {mm}
% V(7)	-	z origin {voxels}
% V(8)	-       y origin {voxels}
% V(9)	-       z origin {voxels}
% W     -       Smoothness {Guassian parameter - voxels}
% FWHM	-       Smoothness {FWHM - mm}
% names -       Sting matrix of parameters in the design matrix
% df    -       degrees of freedom due to error
% Fdf   -       degrees of freedom for the F ratio [Fdf(2) = df]
% TH    -       vector of thresholds used to eliminate extracranial voxels
% FLIP  -       right-left orientation
%	           - 0 left = left  - neurological convention
%	           - 1 left = right - radiological convention
% CONTRAST      - a row vector of contrasts
%
% Image format SPM{Z}
%-----------------------------------------------------------------------
% For each contrast an ANALYZE compatible image format SPM is written to
% a file:  The SPM*.img [and .hdr] are '8 bit' images with the same image
% and voxel sizes as the original data.  Only the positive t values are
% written and the values are scaled by a factor of 16 (i.e. 32 = a Z value
% of 2)
%
% Results matrices in .mat files (at voxels satisfying P{F > f} < UFp)
%-----------------------------------------------------------------------
% XA 	-	adjusted data  		{with grand mean}
% BETA 	-	parameter estimates	{mean corrected}
% XYZ	-	location 		{mm [Talairach]}
% SPMF	-	omnibus F statistic
% SPMt	-	SPM{Z}
%
%__________________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%

% ANALYSIS PROPER
%=======================================================================
global CWD UFp
cd(CWD);

%-Delete files from previous analyses, if they exist
%-----------------------------------------------------------------------
delete XA.mat
delete BETA.mat
delete XYZ.mat
delete SPMF.mat
delete SPMt.mat

%-Location vectors (flipping if the data are radiologically oriented
%-----------------------------------------------------------------------
[y x] = meshgrid([1:V(2,1)],[1:V(1,1)]');
z     = ones(size(x));
x     = (x(:) - ORIGIN(1))*V(4,1);
y     = (y(:) - ORIGIN(2))*V(5,1);
z     =  z(:);

% left right orientation
%-----------------------------------------------------------------------
if FLIP; x = -x; end 			

%-Open image files for SPM{Z}, write headers. SCALE = 1/16, TYPE = 2
%-----------------------------------------------------------------------
for i = 1:size(CONTRAST,1)
	str  = sprintf('SPM%d',i);
	U(i) = fopen([str '.img'],'w');
	spm_hwrite([str '.hdr'],V(1:3,1),V(4:6,1),1/16,2,0,ORIGIN,['spm{Z}']);
end

%-Critical value for F comparison at probability threshold \alpha=UFp
%-----------------------------------------------------------------------
%-NB K, constant partition is only present if H is empty (see spm_spm_ui)
% As we always want to remove mean, a Ko partition is generated for testing
q       = size([K H C B G],1);
r       = rank([K H C B G]);
df      = q - r;
Ko      = ones(q,1);

if ~isempty([H C])
	Fdf = [r - rank([Ko B G]),df];
else
	Fdf = [r - rank(Ko),df];
end
UF      = spm_invFcdf(1 - UFp,Fdf);


%-Initialise variables
%-----------------------------------------------------------------------
D     = zeros(V(1,1)*V(2,1),size(CONTRAST,1)); % dummy matrix for smoothness
sx    = zeros(size(CONTRAST,1),2);             % smoothness estimators {x}
sy    = zeros(size(CONTRAST,1),2);             % smoothness estimators {y}
sz    = zeros(size(CONTRAST,1),2);             % smoothness estimators {z}
TH    = TH*ones(1,V(1,1)*V(2,1));              % global activities
S     = 0;                                     % Volume analyzed
W     = [NaN NaN NaN];


%-Cycle over planes to avoid working memory problems
%-----------------------------------------------------------------------
for i = 1:V(3,1)

  %-Form data matrix for this slice
  %---------------------------------------------------------------------
  X     = zeros(q,V(1,1)*V(2,1));
  for j = 1:q
	tmp    = spm_slice_vol(V(:,j),spm_matrix([0 0 i]),[V(1,1) V(2,1)],0);
	X(j,:) = tmp(:)';
  end

  %-Eliminate background voxels (based on global threshold TH),
  % and eliminate voxels where there are no differences across scans.
  %---------------------------------------------------------------------
  Q = find(all(X > TH) & any(diff(X)));


  if length(Q)
	X     = X(:,Q);
	S     = S + length(Q); 					%-Volume 
	XYZ   = [x(Q) y(Q) z(Q)*(i - ORIGIN(3))*V(6,1)]'; 	%-Locations

	%-if K does not exist remove the grand mean and replace it later.
	% This is simply a device to prevent arbitrary partitoning of
	% the grand mean between the effects modeled in H and B
	%-------------------------------------------------------------------
	if ~size(K,2)
		EX = Ko*mean(X); X = X - EX; end

	%-Estimate parameters and sum of squares due to error
	%-Use pseudo inverse rather than BETA=inv(D'*D)*D'*X for D = DesMtx,
	% to allow for non-unique designs. See matlab help.
	%-------------------------------------------------------------------
	clear BETA
	BETA  = pinv([K H C B G])*X;
	ResSS = sum((X - [K H C B G]*BETA).^2);


	%-Test for overall evidence of effects of interest
	%-Use the likelihood ratio F test - the "extra sum of squares"
	% principle, giving SPM{F} as the F ratio of variances.
	%-------------------------------------------------------------------
	if ~isempty([H C])
		N = sum((X - [Ko B G]*pinv([Ko B G])*X).^2);
	else
		%-Compute F-statistic for confounds & block
		%-(since [H C] is empty there's a K (constant) partition)
		%-Compute adjusted sum of squares, = sum((X - Ko*pinv(Ko)*X).^2)
		N = sum((X - Ko*mean(X)).^2);
	end
	F = (N - ResSS)/Fdf(1)./(ResSS/Fdf(2));


	%-Compute t-statistics for specified compounds of parameters,
	% giving SPM{t}, perform univariate probability transform to z-scores, 
	% giving SPM{Z}.
	%-------------------------------------------------------------------
	T     = zeros(size(CONTRAST,1),size(BETA,2));
	ResMS = ResSS / df;
	for j = 1:size(CONTRAST,1)
		t       = CONTRAST(j,:);
		SE      = sqrt(ResMS*(t*pinv([K H C B G]'*[K H C B G])*t'));
		d       = t*BETA./SE;
		T(j,:)  = spm_t2z(d,df);
		d       = zeros(V(1,1)*V(2,1),1);
		d(Q)    = T(j,:);

        	%-write SPM{Z}
		%-----------------------------------------------------------
        	fwrite(U(j),d.*(d > 0)*16,spm_type(2));

		%-Smoothness estimation
		% Remove the mean of the SPM{Z} prior to smoothness estimation
		%-----------------------------------------------------------
		d       = zeros(V(1,1)*V(2,1),1);
		d(Q)    = T(j,:) - mean(T(j,:));

		%-Compute sums of squares of SPM{Z}, and spatial derivatives
		%-----------------------------------------------------------
		dz      = D(:,j) - d;
		d       = reshape(d,V(1,1),V(2,2));
		[dy dx] = gradient(d);
		Y       = ~d;
		Y       = ~(Y | abs(gradient(Y))); 
		Y       = Y(:);
		sx(j,1) = sx(j,1) + sum( d(Y).^2);
		sx(j,2) = sx(j,2) + sum(dx(Y).^2);
		sy(j,1) = sy(j,1) + sum( d(Y).^2);
		sy(j,2) = sy(j,2) + sum(dy(Y).^2);
		sz(j,1) = sz(j,1) + sum( d(D(:,j) & d(:)).^2);
		sz(j,2) = sz(j,2) + sum(dz(D(:,j) & d(:)).^2);
		D(:,j)  = d(:);
	end % (for j)

	%-Adjustment (remove effects of no interest)
	%-------------------------------------------------------------------
	clear XA
	XA = X - [zeros(size([K H C])) B G]*BETA;
	if ~size(K,2)
		XA = XA + EX; end

	%-Remove voxels with (uncorrected) non-significant F-statistic
	% from further analysis
	%-------------------------------------------------------------------
	Q = find(F > UF);

	%-Cumulate these voxels
	%-------------------------------------------------------------------
	if length(Q)
		spm_append('XA',XA(:,Q));
		spm_append('SPMF',F(Q) );
		spm_append('BETA',BETA(:,Q));
		spm_append('XYZ',XYZ(:,Q));
		if ~isempty(T); spm_append('SPMt',T(:,Q)); end
	end	% (if length(Q))

  else %-No voxels to analyze in this plane

	for j = 1:size(CONTRAST,1)
		fwrite(U(j),zeros(V(1,1)*V(2,1),1),spm_type(2));
	end
  end 			% (conditional on non-zero voxels)
end			% (loop over planes)

%-Smoothness estimates 
%-----------------------------------------------------------------------
if ~isempty(sx)
	W = sqrt([sx(:,1)./sx(:,2) sy(:,1)./sy(:,2) sz(:,1)./sz(:,2)]/2); end
if V(3,1) == 1;   W = W(:,1:2);  end			% 2 dimnesional data
if size(W,1) > 1; W = mean(W); end			% average over contrasts

FWHM  = sqrt(8*log(2))*W.*V(([1:length(W)] + 3),1)'; 	% FWHM in mm


%-Unmap volumes
%-----------------------------------------------------------------------
for i  = 1:q; spm_unmap_vol(V(:,i)); end

%-Save design matrix, and other key variables; S UF CONTRAST W V and FWHM
%-----------------------------------------------------------------------
V      = [V(1:6,1); ORIGIN(:)];
TH     = TH(1,:);
save SPM K H C B G S UF V W FWHM CONTRAST df Fdf TH FLIP names


%-Display and print SPM{F}, Design matrix and textual information
%=======================================================================
figure(3); spm_clf
load XYZ
load SPMF
axes('Position',[-0.05 0.5 0.8 0.4]);
spm_mip(sqrt(SPMF),XYZ,V(1:6))
title(sprintf('SPM{F} p < %f, df: %d,%d',UFp,Fdf),'FontSize',16,'Fontweight','Bold')
text(240,220,sprintf('Search volume: %d voxels',S))
text(240,240,sprintf('Image size: %d %d %d voxels',V(1:3)))
text(240,260,sprintf('Voxel size  %0.1f %0.1f %0.1f mm',V(4:6)))
text(240,280,sprintf('Resolution {FWHM} %0.1f %0.1f %0.1f mm',FWHM))



%-Print out contrasts
%-----------------------------------------------------------------------
axes('Position',[0.1 0.1 0.8 0.4],'XLim',[0,1],'YLim',[0,1]); axis off
line([0 1],[1 1],'LineWidth',3);
line([0 1],[0.92,0.92],'LineWidth',3);
line([0 1],[0.85 0.85],'LineWidth',1);

text(0,0.96,'Results directory:')
text(0.23,0.96,CWD,'FontSize',12,'Fontweight','Bold')
text(0,0.88,'Contrasts','FontSize',12,'Fontweight','Bold')

x0 = 0.25; y0 = 0.82; dx = 0.1; dy = 0.04;
for j = 1 + size(K,2):size([K H C],2)
	text(0,y0 - (j - 1 - size(K,2))*dy,names(j,:),'FontSize',10)
end % (for)
line([0 1],[1,1]*(y0 - size([H C],2)*dy),'LineWidth',1);
for i = 1:size(CONTRAST,1)
	text(x0 + dx*(i - 1),0.88,int2str(i),'FontSize',10)
	for j = 1 + size(K,2):size([K H C],2)
		text(x0 + dx*(i - 1),y0 - (j - 1 - size(K,2))*dy,...
			sprintf('%-6.4g',CONTRAST(i,j)),'FontSize',10)
	end % (for)
end % (for)

spm_print


%-Display, characterize and print SPM{Z}
%-----------------------------------------------------------------------
if ~isempty(CONTRAST)
	load SPMt
	for i = 1:size(CONTRAST,1)
	    spm_projections...
	    (SPMt(i,:),XYZ,Ut,V,W,S,[K H C B G],CONTRAST(i,:),df,1,names);
	    spm_print
	end % (for ...)
end % (if ...)

%-Clear figure 2, the input window
%-----------------------------------------------------------------------
figure(2); clf
set(2,'Name',' ','Pointer','Arrow');
