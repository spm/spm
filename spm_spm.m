function spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,TH,Dnames,Fnames,SIGMA,RT)
% Statistical analysis with the General linear model
% FORMAT spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,TH,Dnames,Fnames,SIGMA,RT);
%
% V   - {12 x q} matrix of identifiers of memory mapped data {q scans}
% H   - {q  x h} condition subpartition of the design matrix {h conditions}
% C   - {q  x c} covariate subpartition of the design matrix {c covariates}
% B   - {q  x n} block     subpartition of the design matrix {n subjects}
% G   - {q  x g} confound  subpartition of the design matrix {g covariates}
%
% CONTRAST - matrix of contrasts, one per row, with p elements.
% ORIGIN   - the voxel correpsonding to [0 0 0] in mm
% TH       - thresholds for each image defining voxels of interest
% Dnames   - string matrix of Dnames for effects in the design matrix
% Fnames   - string matrix of Filenames corresponding to observations
% SIGMA    - Gaussian parameter of K for correlated observations
% RT       - Repeat time for EPI fMRI (generally interscan interval)
%_______________________________________________________________________
%
% spm_spm is the heart of the SPM package and implements the general
% linear model in terms of a design matrix (composed of H C B
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
% analysis; including the partitions of the design matrix [H C B G], the
% number of voxels analyzed (S), image and voxel dimensions [V],
% smoothness estimates of the SPM{Z} [W], threshold for SPM{F}
% [UF] and the contrasts used [CONTRAST].  See below for a complete
% listing.
%
% Output to the results window includes maximum intensity projections
% of the SPM{F}, the design matrix and a series of pages for the SPM{Z}
% (see 'Results' in the help application).
%
% Variables saved in SPM.mat
%-----------------------------------------------------------------------
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
% df    -       degrees of freedom due to error
% Fdf   -       degrees of freedom for the F ratio [Fdf(2) = df]
% TH    -       vector of thresholds used to eliminate extracranial voxels
% SIGMA -       Gaussian parameter of K for correlated observations
% RT    -       Repeat time for EPI fMRI (generally interscan interval)
% Dnames   -    Sting matrix of parameters in the design matrix
% Fnames   -    string matrix of Filenames corresponding to observations
% CONTRAST -    row vectors of contrasts
%
% Results matrices in .mat files (at voxels satisfying P{F > f} < UFp)
%-----------------------------------------------------------------------
% XA 	-	adjusted data  		{with grand mean}
% BETA 	-	parameter estimates	{mean corrected}
% XYZ	-	location 		{mm [Talairach]}
% SPMF	-	SPM{F}
% SPMt	-	SPM{Z}
%
%_______________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%

% ANALYSIS PROPER
%=======================================================================
global UFp


%-Delete files from previous analyses, if they exist
%-----------------------------------------------------------------------
spm_unlink XA.mat BETA.mat XYZ.mat SPMF.mat SPMt.mat


% temporal convolution of the design matrix - dispersion = SIGMA
% the block partition is deliberately omitted here (B is used to associate
% scans and subjects in subsequent routines)
%-----------------------------------------------------------------------
q     = size([H C B G],1);
K     = spm_sptop(SIGMA,q);
H     = K*H;
C     = K*C;
G     = K*G;


%-Location vectors {Xp - voxels, Xq - mm}
%-----------------------------------------------------------------------
x     = 1:V(1); 
y     = 1:V(2);
z     = 1:V(3);
[yp xp] = meshgrid(y,x');
xp    = xp(:)';
yp    = yp(:)';
zp    = zeros(size(yp));
op    = ones(size(yp));
p2mm  = spm_matrix([0 0 0 0 0 0 V(4:6,1)'])*spm_matrix(-ORIGIN);

%-Critical value for F comparison at probability threshold (UFp)
%-----------------------------------------------------------------------
DESMTX = [H C B G];
Fdf    = spm_AnCova([H C],[B G],SIGMA);
UF     = spm_invFcdf(1 - UFp,Fdf);
df     = Fdf(2);


%-Initialise variables
%-----------------------------------------------------------------------
S      = 0;                                     	% Volume analyzed
n_res  = min([q 32]);					% Residuals to be
i_res  = linspace(1,q,n_res);				% used in smoothness
PrevPl = zeros(V(1,1)*V(2,1),n_res);			% estimation
sx_res = 0;              
sy_res = 0;             
sz_res = 0; 
nx_res = 0;
ny_res = 0;
nz_res = 0;            



%-Cycle over planes to avoid working memory problems
%-----------------------------------------------------------------------
spm_progress_bar('Init',100,'AnCova',' ');

for i = 1:V(3)

	%-identify intracranial voxels
	%---------------------------------------------------------------
	Xp    = [xp; yp; (z(i) + zp); op];
	X     = spm_sample_vol(V(:,1),Xp(1,:)',Xp(2,:)',Xp(3,:)',0);
	Q     = find(X > TH(1));
	if length(Q) % proceed


	%-get data
	%---------------------------------------------------------------
	U     = Q;
	X     = zeros(q,length(Q));
	Xp    = Xp(:,Q);
	for j = 1:q
		d      = spm_sample_vol(V(:,j),Xp(1,:)',Xp(2,:)',Xp(3,:)',0);
		U      = U & (d > TH(j));
		X(j,:) = d';
	end
	U     = find(U);
	if length(U); % proceed


	%-volume and locations
	%---------------------------------------------------------------
	Q     = Q(U);
	Xp    = Xp(:,U);
	Y     = X(:,U); clear X
	S     = S + length(Q); 
	XYZ   = p2mm*Xp;
	XYZ   = XYZ([1:3],:);


	%-Remove the grand mean and replace it later
	%---------------------------------------------------------------
	EX    = mean(Y);
	if  q > 512
		for j = 1:size(Y,2); Y(:,j) = Y(:,j) - EX(j); end
	else
		Y     = Y - ones(q,1)*EX;
	end

	%-Convolve
	%---------------------------------------------------------------
	X     = K*Y; clear Y


	%-AnCova; employing pseudoinverse to allow for non-unique designs	
	%---------------------------------------------------------------
	[Fdf F BETA T] = spm_AnCova([H C],[B G],SIGMA,X,CONTRAST);
	Res            = X(i_res,:) - DESMTX(i_res,:)*BETA;	


	%-Remove voxels with (uncorrected) non-significant F-statistic
	%---------------------------------------------------------------
	P     = find(F > UF);
	if length(P) % proceed

	F     = F(P);
	BETA  = BETA(:,P);
	XYZ   = XYZ(:,P);

	%-Adjustment: remove confounds and replace grand mean
	%---------------------------------------------------------------
	XA  = X(:,P) - [zeros(size([H C])) B G]*BETA + ones(q,1)*EX(:,P);


	%-Cumulate remaining voxels
	%---------------------------------------------------------------
	spm_append('XA',XA);
	spm_append('SPMF',F );
	spm_append('BETA',BETA);
	spm_append('XYZ',XYZ);
	if ~isempty(T)
		spm_append('SPMt',spm_t2z(T(:,P),df)); end


	end % proceed P
	clear BETA XA 


	% Normalize residuals
	%---------------------------------------------------------------
	ResSS  = sqrt(sum(Res.^2));
	for  j = 1:size(Res,1)
		Res(j,:) = Res(j,:)./ResSS; end


	%-Compute spatial derivatives
	%---------------------------------------------------------------
	Qx     = zeros(V(1)*V(2),1);
	Qx(Q)  = Q;
	Qx     = reshape(Qx,V(1),V(2));
	Qx     = find(~(~Qx | abs(gradient(~Qx))));
	Qy     = Qx;
	Qz     = zeros(V(1)*V(2),1);
	Qz(Q)  = Q;
	Qz     = find(PrevPl(:,1) & Qz(:));
	nx_res = nx_res + length(Qx);
	ny_res = ny_res + length(Qy);
	nz_res = nz_res + length(Qz);
	for j  = 1:size(Res,1)
		CurrPl      = zeros(V(1)*V(2),1);
		CurrPl(Q)   = Res(j,:);
		dz          = PrevPl(:,j) - CurrPl;
		PrevPl(:,j) = CurrPl;
		[dy dx]     = gradient(reshape(CurrPl,V(1),V(2)));

		d       = dx(Qx); d = d(finite(d));
		sx_res  = sx_res + sum(d.^2);
		d       = dy(Qy); d = d(finite(d));
		sy_res  = sy_res + sum(d.^2);
		d       = dz(Qz); d = d(finite(d));
		sz_res  = sz_res + sum(d.^2);
	end


	end % proceed U
	end % proceed Q

	% progress
	%---------------------------------------------------------------
	spm_progress_bar('Set',i*100/V(3));

end  % (loop over planes)
spm_progress_bar('Clear');


%-Smoothness estimates %-----------------------------------------------------------------------
Lc2z   = spm_lambda(df);
L_res  = [sx_res/nx_res sy_res/ny_res sz_res/nz_res]*(df - 2)/(df - 1);
W      = (2*Lc2z*L_res).^(-1/2);
if V(3) == 1
	W = W(:,[1:2]); end		% 2 dimnesional data

%-Unmap volumes
%-----------------------------------------------------------------------
for i  = 1:q; spm_unmap(V(:,i)); end

%-Save design matrix, and other key variables; S UF CONTRAST W V and df
%-----------------------------------------------------------------------
V      = [V(1:6,1); ORIGIN(:)];
TH     = TH(:,1);
save SPM H C B G S UF V W CONTRAST df Fdf TH Dnames Fnames SIGMA RT


%-Display and print SPM{F}, Design matrix and textual information
%=======================================================================
FWHM   = sqrt(8*log(2))*W.*V(([1:length(W)] + 3))';
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph); spm_clf(Fgraph)
if exist('SPMF.mat')
	load XYZ
	load SPMF
	axes('Position',[-0.05 0.5 0.8 0.4]);
	spm_mip(sqrt(SPMF),XYZ,V)
	str = sprintf('SPM{F} p < %0.2f, df: %0.1f,%0.1f',UFp,Fdf);
	title(str,'FontSize',16)
end
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
text(0.23,0.96,pwd,'FontSize',12,'Fontweight','Bold')
text(0,0.88,'Contrasts','FontSize',12,'Fontweight','Bold')

x0 = 0.25; y0 = 0.82; dx = 0.8/size(CONTRAST,1); dy = 0.04;
for j = 1:size([H C],2)
	text(0,y0 - (j - 1)*dy,Dnames(j,:),'FontSize',10)
end

line([0 1],[1,1]*(y0 - size([H C],2)*dy),'LineWidth',1);
for i = 1:size(CONTRAST,1)
	text(x0 + dx*(i - 1),0.88,int2str(i),'FontSize',10)
	for j = 1:size([H C],2)
		str = sprintf('%-6.3g',CONTRAST(i,j));
		text(x0 + dx*(i - 1),y0 - (j - 1)*dy,str,'FontSize',10)
	end
end

spm_print


%-Display, characterize and print SPM{Z}
%-----------------------------------------------------------------------
if exist('SPMt.mat')
	load SPMt
	U     = spm_invNcdf(1 - 0.01);
	[P,EN,Em,En] = spm_P(1,W,U,0,S);
	K     = round(En);
	for i = 1:size(CONTRAST,1)
	    spm_projections(SPMt(i,:),XYZ,U,K,V,W,S,DESMTX,CONTRAST(i,:),df);
	    spm_print
	end
end

