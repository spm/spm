function spm_fmri_spm(V,H,C,B,G,CONTRAST,ORIGIN,GX,RT,SIGMA)
% Statistical analysis with the General linear model
% FORMAT spm_fmri_spm(V,H,C,B,G,CONTRAST,ORIGIN,GX,RT,SIGMA);
%
% V   - {12 x q} matrix of identifiers of memory mapped data {q scans}
%
% H   - {q  x h} Factor    subpartition of the design matrix {h levels}
% C   - {q  x c} covariate subpartition of the design matrix {c covariates}
% H   - {q  x h} subject   subpartition of the design matrix {b levels}
% G   - {q  x g} confound  subpartition of the design matrix {g covariates}
%
% CONTRAST  - matrix of (c) contrasts {one per row with h + c + n + g elements}
% ORIGIN    - the voxel correpsonding to [0 0 0] in mm
% GX        - Global activities
% RT        - Repeat time
% SIGMA     - temporal smoothing if SIGMA == 1
%____________________________________________________________________________
%
% spm_fmri_spm is the heart of the SPM package and implemets the general
% linear model in terms of a design matrix (composed of H C B and G)
% and the data (V).  Significant compounds of the estimated parameters are
% assessed with a quotient that has the t distribution under the null
% hyypothesis.  The resulting SPM{t} is transformed to the Unit Gaussian
% distribution [SPM{Z}] and characterized by further analysis using the theory 
% of Gaussian Fields (see spm_projections.m for more details)
%
% H, C, B and G are partitions of the design matrix corresponding to covariates
% of interest (e.g. a reference waveform) and factors/covariates of no interest
% (confounds).  The difference between this routine and more conventional
% implementations of the general linear model is that the data (and
% design matrix) are smoothed over observations (time).  The resulting
% error terms have a non-independent, stationary covariance structure that
% is reflected in the 'effective degrees of freedom', that are less than
% the number of scans.  The temporal smoothing uses a Gaussian kernel with
% a standard deviation of sqrt(8) seconds.
%
% The outputs of this routine are a series of .mat files containing
% paramter estimates, adjusted values, SPM{Z} etc that are written to
% CWD (see spm_defaults.m).  IMPORTANT: Existing results are overwritten
% without prompting
%
% Voxels are retained for further analysis if the F ratio for that voxel is
% significant (p < 0.05 uncorrected) and all the voxels have a reasonably 
% high level of activity (0.8 of the global activity).
%
%  SPMF.mat contains a 1 x N vector of F values reflecting the omnibus 
% significance of effects [of interest] at each of the N 'significant' voxels. 
% 'Significance' is defined by the p-value of the F threshold (p < 0.05).
%  XYZ.mat contains a 3 x N matrix of the x,y and z location of the voxels in 
%  SPMF in mm (usually referring the the standard anatomical space (Talairach
% and Tournoux 1988)} (0,0,0) corresponds to the centre of the voxel specified 
% by ORIGIN in the *.hdr of the original and related data.
%   BETA.mat contains a p x N matrix of the p parameter estimates at each of
% the N voxels.  These parameters include all effects specified by the design 
% matrix.
%   XA.mat contains a q x N matrix of adjusted activity values having removed 
% the effects of no interest at each of the N voxels for all q scans.
%   SPMt.mat contains a c x N matrix of the c SPM{Z} defined by the c 
% contrasts supplied for all N voxels at locations XYZ.
%   SPM.mat  contains a collection of matrices that pertain to the analysis; 
% including the partitions of the deign matrix [H C B G], the number of voxels 
% analyzed (S), image and voxel dimensions [V], smoothness estimates of 
% the SPM{Z} [W FWHM], threshold for SPM{F} [UF] and the contrasts [CONTRAST]
%
% Output to the results window includes maximum intensity projections of the 
% SPM{F}, the design matrix and a series of pages for the SPM{Z} (see 
% 'Projections' in the help application).
%
%
% Variables saved in SPM.mat
%---------------------------------------------------------------------------
% H	-	factor    partition of design matrix (usually empty)
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
% V(8)	-	y origin {voxels}
% V(9)	-	z origin {voxels}
% W	-	Smoothness {Guassian parameter - voxels}
% FWHM	-	Smoothness {FWHM - mm}
% GX    -       whole brain means for each scan
% df    -       effective degrees of freedom
% Fdf   -       degrees of freedom for the F ratio
% names -       names of effects in design matrix (usually empty)
% RT    -       Interscan Interval
% FLIP  -       Left-right orientation
%	           - 0 left = left  - neurological convention
%	           - 1 left = right - radiological convention (default)
% CONTRAST      - row matrix of contrasts or compound coefficients
%---------------------------------------------------------------------------
%
% Image format SPM{Z}
%
% For each contrast an ANALYZE compatible image format SPM is written to
% a file:  The SPM*.img [and .hdr] are 'float' images with the same image
% and voxel sizes as the original data.
%---------------------------------------------------------------------------
%
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% - in press
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage submitted
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% in press
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage - in press
%
%__________________________________________________________________________
% %W% Karl Friston %E%
 


% results matrices in .mat files
%---------------------------------------------------------------------------
% XA 	-	adjusted data  		{with at voxels F: p < 0.05}
% BETA 	-	parameter estimates	{mean corrected}
% XYZ	-	location 		{mm [Talairach]}
% SPMF	-	omnibus F statistic
% SPMt	-	SPM{Z}



% ANALYSIS PROPER
%===========================================================================
global CWD

% radiological convention (image left = subject right) assumed for T2* data
%---------------------------------------------------------------------------
FLIP  = (V(3) ~= 1);				

% degrees of freedom (specifed by the sizes of the design matrix partitions)
%---------------------------------------------------------------------------
q     = size([H C B G],1);				% total + mean
c     = rank([H C]);					% covariates
g     = rank([B G]);					% confounds


% temporal convolution kernel - dispersion = 8 seconds [delay = sqrt(8)]
%---------------------------------------------------------------------------
K     = eye(q);
if SIGMA
	sigma = sqrt(8)/RT;
	K     = toeplitz( exp(-[0:(q - 1)].^2/(2*sigma^2)) );
	K     = K./(sum(K)'*ones(1,q));
end
H     = K*H;
C     = K*C;
G     = K*G;
% the block partition is deliberately omitted here (B is used to associate
% scans and subjects in subsequent routines)

% the [effective] degrees of freedom -
%---------------------------------------------------------------------------
KK    = K*K';
IG    = pinv([H C B G]);
R     = eye(q) - [H C B G]*IG;
J     = IG*KK*IG';
ER    = trace(R*KK);
df    = round(ER^2/trace(R*KK*R*KK));			% residuals [error]


% Point list of locations {Xp} voxels
%---------------------------------------------------------------------------
x     = 1:V(1,1); 
y     = 1:V(2,1);
z     = 1:V(3,1);

xp    = x'*ones(size(y));
yp    = ones(size(x'))*y;
zp    = zeros(size(yp));

N     = length(x)*length(y);
d     = ones(1,N);
Xp    = zeros(4, length(x)*length(y)*length(z));
Xq    = zeros(3, length(x)*length(y)*length(z));
for i = 1:length(z)
	j  = (i - 1)*N + [1:N];
	Xp(:,j) = [xp(:)'; yp(:)'; (z(i) + zp(:)'); d];
end

% and {Xq} in mm
%---------------------------------------------------------------------------
Xq      = spm_matrix([0 0 0 0 0 0 V(4:6,1)'])*spm_matrix(-ORIGIN)*Xp;
Xq(4,:) = [];


% left right orientation
%---------------------------------------------------------------------------
if FLIP; Xq(1,:) = -Xq(1,:); end 			

% variables saved (at voxels satisfying P{F > x} < 0.05)
%---------------------------------------------------------------------------
eval(['cd ' CWD]);
delete XA.mat
delete BETA.mat
delete XYZ.mat
delete SPMF.mat
delete SPMt.mat

% open image files for SPM{Z} and write the headers. SCALE = 1, TYPE = 16
%---------------------------------------------------------------------------
for i = 1:size(CONTRAST,1)
	d    = sprintf('SPM%0.0f',i);
	U(i) = fopen([d '.img'],'w');
	spm_hwrite([d '.hdr'],V(1:3,1),V(4:6,1),1,16,0,ORIGIN,['spm{Z}']);
end

% set thresholds
%----------------------------------------------------------------------------
if c
   UF = spm_invFcdf(1 - 0.05,[c,df]);		% threhold for SPM{F}
else
   UF = spm_invFcdf(1 - 0.05,[g,df]);
end
Ut    = spm_invNcdf(1 - 0.01);			% threhold for SPM{t}/SPM{Z}
S     = 0;					% Volume analyzed
dI    = 256;					% number of voxels per cycle


% cycle over 'sets' of dI voxels to avoid working memory problems
%---------------------------------------------------------------------------
for i   = 0:dI:length(Xp)

  I     = [1:dI] + i;
  I     = I(I <= length(Xp));
  X     = zeros(q,length(I));				% data matrix
  for j = 1:q
      d = spm_sample_vol(V(:,j),Xp(1,I)',Xp(2,I)',Xp(3,I)',0);
      X(j,:) = d(:)';
  end
  XYZ   = Xq(:,I);					% locatinns


  % eliminate background voxels (less than 0.8 of the global activity GX)
  %-------------------------------------------------------------------------
  Q     = find(all(X > 0.8*GX*ones(1,size(X,2))) & ~all(~diff(X)));


  if length(Q)
    X     = X(:,Q);
    XYZ   = XYZ(:,Q);						% locatinns
    S     = S + length(Q);					% volume so far

    % convolution
    %-----------------------------------------------------------------------
    X     = K*X;

    % estimate parameters and sum of squares due to error
    %-----------------------------------------------------------------------
    BETA  = IG*X;
    R     = sum((X - [H C B G]*BETA).^2);


    % test for effects of interest SPM{F} with he F ratio of variances
    %-----------------------------------------------------------------------
    if isempty([B G]) | isempty([H C])
	N = sum(X.^2);
	n = c + g;
    else
	N = sum((X - [B G]*pinv([B G])*X).^2);
	n = c;
    end
    F     = (N - R)/n./(R/df);

    % test of linear compound of parameters SPM{t} and transform to SPM{Z}
    %-----------------------------------------------------------------------
    T     = zeros(size(CONTRAST,1),size(BETA,2));
    E     = R/ER;
    for j = 1:size(CONTRAST,1)
	t       = CONTRAST(j,:);
	RMS     = sqrt(E*(t*J*t'));
	d       = t*BETA./RMS;

	% transform and cumulate distribution
        %-------------------------------------------------------------------
	Z       = spm_t2z(d,df);
	T(j,:)  = Z;
	d       = zeros(length(I),1);
	d(Q)    = T(j,:);

        % write SPM{Z}
        %-------------------------------------------------------------------
        fwrite(U(j),d,spm_type(16));

    end

    % adjustment (remove effects of no interest)
    %-----------------------------------------------------------------------
    d      = ([1:size([B G],2)]  + size([H C],2));
    if isempty(d)
	XA = X;
    else
	XA = X - [B G]*BETA(d,:);
    end

    % remove non-significant voxels from further analysis
    %-----------------------------------------------------------------------
    Q      = find(F > UF);

    % cumulate these voxels
    %-----------------------------------------------------------------------
    if length(Q)
    	spm_append('XA'  ,XA(:,Q));
    	spm_append('SPMF',F(Q) );
    	spm_append('BETA',BETA(:,Q));
    	spm_append('XYZ' ,XYZ(:,Q));
    	if ~isempty(T ); spm_append('SPMt',T(:,Q)); end

     end				% end conditional on P(F > x) < 0.05
  else
    	for j = 1:size(CONTRAST,1)
        	fwrite(U(j),zeros(length(I),1),spm_type(16));
	end
  end					% end conditional on non-zero voxels
end					% end cycle over planes

% smoothness estimates %---------------------------------------------------------------------------
W     = [];
FWHM  = [];
for i = 1:size(CONTRAST,1)
	[d fwhm] = spm_W(spm_map(sprintf('SPM%0.0f.img',i)));
	W        = [W; d];
	FWHM     = [FWHM; fwhm];
end
if size(W,1) > 1
	W        = mean(W);
	FWHM     = mean(FWHM);
elseif size(W,1) < 1
	W        = [NaN NaN NaN];
	FWHM     = W;
end


% unmap volumes and save design matrix S UF CONTRAST W V and FWHM
%---------------------------------------------------------------------------
for i  = 1:q; spm_unmap_vol(V(:,i)); end

% save design matrix and other key variables
%---------------------------------------------------------------------------
V      = [V(1:6,1); ORIGIN(:)];
Fdf    = [rank([H C]) df];
names  = [];

save SPM H B C G S UF V W FWHM CONTRAST GX sigma df Fdf RT names FLIP


% display and print SPM{F}, Design matrix and textual information
%---------------------------------------------------------------------------
figure(3); spm_clf

load XYZ
load SPMF
axes('Position',[0.1 0.5 0.8 0.4]);
spm_mip(sqrt(SPMF),XYZ,V(1:6))
title('SPM{F} (p < 0.05)','FontSize',16,'Fontweight','Bold')

axes('Position',[0.1 0.1 0.36 0.36]);
imagesc(spm_DesMtxSca([H B C G]))
axis image
title ('Design matrix','Fontweight','Bold')
xlabel 'effect'; ylabel 'scan'


axes('Position',[0.46 0.1 0.42 0.3]); axis off
text(0,1.05,'Results directory:')
text(0,.95,CWD,'FontSize',16,'Fontweight','Bold')
text(0,.8,sprintf('Search volume %0.0f voxels',S))
text(0,.7,sprintf('Resolution {FWHM} %-4.1f %-4.1f %-4.1f mm',FWHM))
text(0,.6,sprintf('Image size  %-4.0f %-4.0f %-4.0f voxels',V(1:3)))
text(0,.5,sprintf('Voxel size  %-4.1f %-4.1f %-4.1f mm',V(4:6)))
text(0,.35,'Degrees of freedom','Fontsize',10)
text(0,.25,sprintf('Covariates  %0.0f ',rank([H C])),'Fontsize',10)
text(0,.20,sprintf('Confounds  %0.0f ',rank([B G])),'Fontsize',10)
text(0,.15,sprintf('Residuals  %0.0f ',df),'Fontsize',10)
line([0 1],[1.2 1.2],'LineWidth',3);
line([0 1],[.85 .85],'LineWidth',3);
line([0 1],[0.3 0.3],'LineWidth',1);
line([0 1],[0.0 0.0],'LineWidth',1);

spm_print



% display, characterize and print SPM{Z}
%---------------------------------------------------------------------------
if ~isempty(CONTRAST)
	load SPMt
	for i = 1:size(CONTRAST,1)
	    spm_projections(SPMt(i,:),XYZ,Ut,V,W,S,[H C B G],CONTRAST(i,:),df);
	    spm_print
	end
end

% clear figure 2
%---------------------------------------------------------------------------
figure(2); clf
set(2,'Name',' ','Pointer','Arrow');
