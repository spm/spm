function spm_spm(VY,xX,xM,F_iX0,varargin)
% Estimation of the General Linear Model
% FORMAT spm_spm(VY,xX,xM,F_iX0,<extra parameters for SPM.mat>)
%
% VY    - nScan x nVar struct array of mapped image volumes
%         Images must have the same orientation, voxel size and data type
%       - Any scaling should have already been applied via the image handle
%         scalefactors.
%
% xX    - Structure containing design matrix information
%       - Required fields are:
%         xX.X      - Design matrix (raw, not temporally smoothed)
%         xX.Xnames - cellstr of parameter names corresponding to columns
%                     of design matrix
%       - Optional fields are:
%         xX.K      - Sparse temporal smoothing matrix (see spm_sptop)
% 		    - Or cell of session-specific structures (see spm_filter)
%                   - Design & data are smoothed using K
%                     (model is K*Y = K*X*beta + K*e)
%                   - Note that K should not smooth across block boundaries
%                   - defaults to speye(size(xX.X,1))
%         xX.xVi    - structure describing intrinsic temporal auto-correlation
%                   - required fields are
%                     xX.xVi.Vi - intrinsic temporal auto-correlation matrix
%                               - defaults to speye(size(xX.X,1))
%                               - specifying a cell array of contraints
%                                 invokes spm_reml to estimate hyperparameters
%                                 assuming Vi is constant over voxels.
%
% xM    - Structure containing masking information, or a simple column vector
%         of thresholds corresponding to the images in VY.
%       - If a structure, the required fields are:
%         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1=>implicit zero/NaN mask)
%         xM.VM - struct array of mapped explicit mask image volumes
% 		- (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or data
%                 type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%       - Note that voxels with constant data (i.e. the same value across
%         scans) are also automatically masked out.
%
% F_iX0 - Indicies of design matrix columns to form the reduced design
%         matrix. The ensuing F-test (implemented through F-contrasts)
%         is labelled "effects of interest", and is to filter the saved
%         data pointlist (Y.mad).
%       - If F_iX0 is a structure array
%	  F_iX0.iX0  - indices as above
%	  F_iX0.name - name of F-contrast
%       - See Christensen for details of F-contrasts
%
% varargin - Additional arguments are passed through to the SPM.mat file
%            with their original names, permitting the caller (UI) function
%            to pass information to the results section.
%
% In addition, the global default UFp is used to set the critical
% F-threshold for pointlist saving of raw data. UFp is an upper tail
% probability threshold (UFp in [0,1]) UFp==1 => no F-filtering (all
% in-mask voxels have their data written to the pointlist Y.mad file,
% which takes a while!) UFp==0 => no pointlist data written. 0<UFp<1 =>
% data written for voxels with F-statistic significant (uncorrected) at
% UFp level.
%
%_______________________________________________________________________
%
% spm_spm is the heart of the SPM package. Given image files and a
% General Linear Model, it estimates the model parameters, residual
% variance, and smoothness of the standardised residual fields, writing
% these out to disk in the current working directory for later
% interrogation in the results section. (NB: Existing analyses in the
% current working directory are overwritten.)
%
% The model is expressed via the design matrix (xX.X). The basic model
% at each voxel is of the form is Y = X*beta * e, for data Y, design
% matrix X, (unknown) parameters B, and residual errors e. The errors
% are assummed to have a normal distribution ~ N(0,sigma^2*Vi) where Vi can
% be specified in xX.xVi.  If the intrinsic correlations Vi are not
% specified Vi defaults to the identity matrix (i.e. i.i.d assumptions).
% The parameters are estimated by ordinary least squares, using matrix
% methods, and the residual variance estimated.  
%
% Note that only a single variance component is considered. This is
% *always* a fixed effects model. Random effects models can be effected
% using summary statistics, summarising the data such that the fixed
% effects in the model on the summary data are those effects of
% interest (i.e. the population effects). In this case the residual
% variance for the model on the summary data contains the appropriate
% variance components in the appropriate ratio. See spm_RandFX.man for
% further details and below.
%
% Under the additional assumption that the standardised residual images
% are non-stationary standard Gaussian random fields, results from
% Random field theory can be applied to estimate the significance
% statistic images (SPM's) accounting for the multiplicity of testing
% at all voxels in the volume of interest simultaneously. The required
% parameters are the volume, and Lambda, the variance covariance matrix
% of partial derivatives of the standardised error fields. 
% spm_est_smoothness estimates the variances of the partial derivatives 
% in the axis directions (the diagonal of Lambda). The covariances (off
% diagonal elements of Lambda) are assumed to be zero.
% 
%                           ----------------
%
% In some instances the i.i.d. assumptions about errors do not hold.  For
% example, with serially correlated (fMRI) data or correlations among the
% levels of a factor in repeated measures designs.  This non-sphericity
% can be specified in terms of constraints (xX.xVi.Vi) Covariance
% components and correlations will be estimated with an ReML (restricted
% maximum likelihood) algorithm using these contraints.  This estimation
% assumes the same correlation structure for each voxel.  The ReML
% estimates are then used to correct for non-sphericity during inference
% by adjusting the statistics and degrees of freedom appropriately.
% Because spm_spm uses an OLS estimator (as opposed to a Gauss-Markov
% estimator) the parameter estimates will not be affected. Note that the
% non-sphereicity correction can give non-integer degrees of freedom.
%
% The volume analsed is the intersection of the threshold masks,
% explicit masks and implicit masks. See spm_spm_ui for further details
% on masking options.
%
% 
%-----------------------------------------------------------------------
%
% The output of spm_spm takes the form of an SPM.mat file of analysis
% parameters, and 'float' Analyze images of the parameter and variance
% estimates. An 8bit zero-one mask image indicating the volume assessed
% is also written out, with zero indicating voxels outside tha analysed
% volume.
%
% In addition, voxels with F-test for "effects of interest" significant
% at UFp have their raw data written out in compressed pointlist format
% to Y.mad (see spm_append for information on the MAD file format).
% Yidx.mat contains a single vector Yidx of length equal to the column
% dimension of Y.mad, whose entries  indicate the XYZ coordinates for
% the correponding columns of Y.mad by indexing the columns of XYZ
% saved in SPM.mat. The Y.mad & Yidx.mat files are used for plotting,
% and are written to allow the results to be self sufficient (i.e. the
% input data can be archived). If a voxel is not represented in the
% XYZ.mat & Y.mad files, then you can't plot it!  (Unless you do
% manually append data to the Y.mad and Yidx.mat files yourself.) (Note
% that the parameter & variance images are written for all voxels
% within the supplied mask, an advance on SPM94/5/6 which only wrote
% results information for voxels surviving the F-threshold.  The
% smoothness and volume analysed are computed for the entire analysis
% volume, as always.)
%
%                           ----------------
%
% SPM.mat                                         - analysis parameters
% The SPM.mat file contains the following variables:
%     SPMid     - string identifying SPM & program version
%     VY        - as input
%
%     xX        - design matrix structure
%               - all the fields of the input xX are retained
%               - the following additional fields are added:
%     xX.V      - V matrix (xX.K*Vi*xX.K')
%     xX.xKXs   - space structure for K*X (xX.K*xX.X), the temporally smoothed
%                 design matrix
%               - given as spm_sp('Set',xX.K*xX.X) - see spm_sp for details
%                 (so xX.xKXs.X is K*X)
%     xX.pKX    - pseudoinverse of K*X (xX.K*xX.X), computed by spm_sp
%                 (this was called PH in spm_AnCova)
%     xX.pKXV   - pinv(K*X)*V  (this was called PHV in spm_AnCova)
%     xX.Bcov   - pinv(K*X)*V*pinv(K*X)' - variance-covariance matrix
%                 of parameter estimates (when multiplied by ResMS)
%                 (NB: BCOV in SPM96 was xX.Bcov/xX.trRV, so that BCOV    )
%                 (multiplied by ResSS was the variance-covariance matrix )
%                 (of the parameter estimates. (ResSS/xX.trRV = ResMS)    )
%     xX.trRV   - trace of R*V, computed efficiently by spm_SpUtil
%     xX.trRVRV - trace of RVRV (this was called trRV2 in spm_AnCova)
%     xX.erdf   - effective residual degrees of freedom (trRV^2/trRVRV)
%     xX.nKX    - design matrix (xX.xKXs.X) scaled for display
%                 (see spm_DesMtx('sca',... for details)
% 
%     xM        - as input (but always in the structure format)
%
%     M         - 4x4 voxel->mm transformation matrix
%     DIM       - image dimensions {voxels} - column vector
%     Vbeta     - cellstr of beta image filenames (relative)
%     VResMS    - relative filename of ResMS image
%     XYZ       - 3xS (S below) vector of in-mask voxel coordinates
%     F_iX0     - structure array of default F-contrasts
%     F_iX0.iX0 - as input
%     F_iX0.name- as input
%     UFp       - critical p-value for F-thresholding
%     UF        - F-threshold value
%     S         - Lebesgue measure or volume (in voxels)
%     R         - vector of resel counts     (in resels)
%     FWHM      - Smoothness (of component fields - FWHM, in voxels)
%                 
%
% ...plus any additional arguments passed to spm_spm for saving in the
%    SPM.mat file
%
%                           ----------------
%
% xCon.mat                                        - contrast definitions
% The xCon.mat file contains a single structure defining the default F-contrast
% (The results section can be used to define additional contrasts.)
%     xCon       - Contrast structure (created by spm_FcUtil.m)
%     xCon.name  - Name of contrast
%     xCon.STAT  - 'F', 'T' or 'X' - for F/T-contrast ('X' for multivariate)
%     xCon.c     - (F) Contrast weights
%     xCon.X0    - Reduced design matrix (spans design space under Ho)
%                  It is in the form of a matrix (spm99b) or the
%                  coordinates of this matrix in the orthogonal basis
%                  of xX.X defined in spm_sp. 
%     xCon.iX0   - Indicates how contrast was specified:
%                  If by columns for reduced design matrix then iX0 contains the
%                  column indices. Otherwise, it's a string containing the
%                  spm_FcUtil 'Set' action: Usuall one of {'c','c+','X0'}
%                  (Usually this is the input argument F_iX0.)
%     xCon.X1o   - Remaining design space (orthogonal to X0).
%                  It is in the form of a matrix (spm99b) or the
%                  coordinates of this matrix in the orthogonal basis
%                  of xX.X defined in spm_sp.
%     xCon.eidf  - Effective interest degrees of freedom (numerator df)
%     xCon.Vcon  - ...for handle of contrast/ESS image (empty at this stage)
%     xCon.Vspm  - ...for handle of SPM image (empty at this stage)
%
%                           ----------------
%
% mask.{img,hdr}                                   - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
%
%                           ----------------
%
% beta_????.{img,hdr}                                 - parameter images
% These are 16-bit (float) images of the parameter estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.img) are given
% value NaN.
%
%                           ----------------
%
% ResMS.{img,hdr}                    - estimated residual variance image
% This is a 32-bit (double) image of the residual variance estimate.
% Voxels outside the analysis mask are given value NaN.
%
%                           ----------------
%
% RVP.{img,hdr}                      - estimated resels per voxel image
% This is a 32-bit (double) image of the RESELs per voxel estimate.
% Voxels outside the analysis mask are given value 0.  These images
% reflect the nonstationary aspects the spatial autocorrelations.
%
%                           ----------------
%
% Yidx.mat [optional]     - containing 1xn vector Yidx of column indices
% If raw data is being written out in compressed pointlist format to a
% Y.mad file (i.e. UFp > 0), then Yidx.mat contains Yidx, a 1xn matrix of
% column indices. XYZ(:,Yidx) are then the voxel coordinates of the
% data in successive columns of Y.mad.
%
%                           ----------------
%
% Y.mad [optional]                    - compressed pointlist of raw data
% The Y.mad file contains the raw (unsmoothed) data at voxels with an
% F-statistic significant at UFp. Essentially, columns contain the data
% at a single voxel, rows correspond to the input images (VY). The
% corresponding columns of XYZ specify the voxel locations.
%
% The MAD format is a compressed format: Each voxel's data is
% individually scaled and shifted to fit into a given range (usually
% 0-2^8 for 8-bit MAD files), and saved using relatively few bits. The
% scale and location parameters are saved as floats. This means it can
% take a while to write a lot of data. However, when plotting, usually
% only data at a single voxel is required, so reading os OK, This
% presents a good compromise between data integrity and file size.
%
% See spm_append & spm_extract for MAD file manipulation.
%
%-----------------------------------------------------------------------
%
% References:
%
% Christensen R (1996) Plane Answers to Complex Questions
%       Springer Verlag
%
% Friston KJ, Holmes AP, Worsley KJ, Poline JP, Frith CD, Frackowiak RSJ (1995)
% ``Statistical Parametric Maps in Functional Imaging:
%   A General Linear Approach''
%       Human Brain Mapping 2:189-210
%
% Worsley KJ, Friston KJ (1995)
% ``Analysis of fMRI Time-Series Revisited - Again''
%       Neuroimage 2:173-181
%
%-----------------------------------------------------------------------
%
% For those interested, the analysis proceeds a "block" at a time,
% The block size conforms to maxMem that can be set as a global variable
% MAXMEM (in bytes) [default = 2^20]
%
%_______________________________________________________________________
% %W% Andrew Holmes, Jean-Baptiste Poline, Karl Friston %E%
SCCSid   = '%I%';


%-Say hello
%-----------------------------------------------------------------------
SPMid    = spm('FnBanner',mfilename,SCCSid);
Finter   = spm('FigName','Stats: estimation...'); spm('Pointer','Watch')


%-maxMem is the maximum amount of data that will be processed at a time
%-----------------------------------------------------------------------
global MAXMEM 
if length(MAXMEM)	
	maxMem = MAXMEM;	
else
	maxMem = 2^20;	%-Max data block size, in bytes
end
Def_UFp  = 0.001;	%-Default F-threshold for Y.mad pointlist filtering
maxRes   = 64;		%-Maximum # res images for smoothness estimation


%-Condition arguments
%-----------------------------------------------------------------------
if nargin < 2, error('Insufficient arguments'), end
if nargin < 3, xM = zeros(size(X,1),1); end
if nargin < 4, c  = []; end


%-Check required fields of xX structure exist - default optional fields
%-----------------------------------------------------------------------
for tmp = {'X','Xnames'}
	if ~isfield(xX,tmp)
     		error(sprintf('xX doesn''t contain ''%s'' field',tmp{:}))
	end
end
if ~isfield(xX,'K')
	xX.K  = speye(size(xX.X,1));
end
if ~isfield(xX,'xVi')
	xX.xVi = struct(	'Vi',	speye(size(xX.X,1)),...
				'form',	'none'); 
end


%-If xM is not a structure then assumme it's a vector of thresholds
%-----------------------------------------------------------------------
if ~isstruct(xM), xM = struct(	'T',	[],...
				'TH',	xM,...
				'I',	0,...
				'VM',	{[]},...
				'xs',	struct('Masking','analysis threshold'));
end

%-Delete files from previous analyses
%-----------------------------------------------------------------------
if exist(fullfile('.','SPM.mat'),'file')==2
	spm('alert!',{...
		'Current directory contains some SPMstats results files!',...
    		['        (pwd = ',pwd,')'],...
		'Existing results are being overwritten!'},...
		mfilename,1);
	warning(sprintf('Overwriting existing results\n\t (pwd = %s) ',pwd))
	drawnow
end

files = {	'SPM.mat','Yidx.mat','Y.mad','xCon.mat',...
		'mask.???','ResMS.???','RVP.???',...
		'beta_????.???','con_????.???','ResI_????.???',...
		'ess_????.???', 'spm?_????.???'};

for i=1:length(files)
	if any(files{i} == '*'|files{i} == '?' )
		[tmp,null] = spm_list_files(pwd,files{i});
		for i=1:size(tmp,1)
			spm_unlink(deblank(tmp(i,:)))
		end
	else
		spm_unlink(files{i})
	end
end

%-Check & note Y images dimensions, orientation & voxel size
%-----------------------------------------------------------------------
if any(any(diff(cat(1,VY(:).dim),1,1),1) & [1,1,1,0]) 
	error('images do not all have the same dimensions')
end
if any(any(any(diff(cat(3,VY(:).mat),1,3),3)))
	error('images do not all have same orientation & voxel size')
end

M      = VY(1,1).mat;
DIM    = VY(1,1).dim(1:3)';
N      = 3 - sum(DIM == 1);


%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Initialise design space
%=======================================================================
fprintf('%-40s: %30s','Initialising design space','...computing')    %-#

%-Construct design parmameters, and store in design structure xX
% Take care to apply temporal confounds - KX stored as xX.xKX.X
%-Note that Vi may not be known exactly at this point, if it is to be
% estimated. Parameters dependent on Vi are committed to xX at the end.
%-Note that the default F-contrast (used to identify "interesting" voxels
% to save raw data for) computation requires Vi. Thus, if Vi is to be
% estimated, any F-threshold will only have upper tail probability UFp. 
%-----------------------------------------------------------------------
% V            - Autocorrelation matrix K*Vi*K'
% xX.xKXs      - Design space structure of KX
% xX.pKX       - Pseudoinverse of KX
% trRV,trRVRV  - Variance expectations
% erdf         - Effective residual d.f.
%-----------------------------------------------------------------------

%-Get provisional correlation structure (V)
%-----------------------------------------------------------------------
[nScan nBeta] = size(xX.X);
[nScan nVar]  = size(VY);
nSres         = min(nScan,maxRes);

if iscell(xX.xVi.Vi)
	Vi    = speye(nScan);
else
	Vi    = xX.xVi.Vi;
end
KVi           = spm_filter(xX.K, Vi);
V             = spm_filter(xX.K,KVi');

%-Parameter projection matrix and traces
%-----------------------------------------------------------------------
xX.xKXs       = spm_sp('Set',spm_filter(xX.K, xX.X));
xX.pKX        = spm_sp('x-',xX.xKXs);
[trRV trRVRV] = spm_SpUtil('trRV',xX.xKXs,V);
erdf          = trRV^2/trRVRV;


%-Check estimability
%-----------------------------------------------------------------------
if  erdf < 0
    error(sprintf('This design is completely unestimable! (df=%-.2g)',erdf))
elseif erdf == 0
    error('This design has no residuals! (df=0)')
elseif erdf <  4
    warning(sprintf('Very low degrees of freedom (df=%-.2g)',erdf))
end


%-Default F-contrasts (in contrast structure) & Y.mad pointlist filtering
%=======================================================================
fprintf('%s%30s',sprintf('\b')*ones(1,30),'...F-contrast')           %-#
UFp   = spm('GetGlobal','UFp'); if isempty(UFp), UFp = Def_UFp; end

if isempty(F_iX0)
	F_iX0 = struct(	'iX0',		[],...
			'name',		'all effects');
elseif ~isstruct(F_iX0)
	F_iX0 = struct(	'iX0',		F_iX0,...
			'name',		'effects of interest');
end

%-Create Contrast structure array
%-----------------------------------------------------------------------
xCon  = spm_FcUtil('Set',F_iX0(1).name,'F','iX0',F_iX0(1).iX0,xX.xKXs);
for i = 2:length(F_iX0)
	xcon = spm_FcUtil('Set',F_iX0(i).name,'F','iX0',F_iX0(i).iX0,xX.xKXs);
	xCon = [xCon xcon];
end

%-Parameters for saving in Y.mad (based on first F-contrast)
%-----------------------------------------------------------------------
[trMV trMVMV] = spm_SpUtil('trMV',spm_FcUtil('X1o',xCon(1),xX.xKXs),V);
eidf          = trMV^2/trMVMV;
h             = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);


%-Compute UF, the F-threshold
%-----------------------------------------------------------------------
if UFp > 0 & UFp < 1				%-F-filter for Y.mad file
	UF   = spm_invFcdf(1 - UFp,[eidf,erdf]);
	Yidx = [];
elseif UFp == 1					%-No filtering - save all data
	UF   = -Inf;
	Yidx = [];
elseif UFp == 0					%-Write no Y.mad data
	UF   =  Inf;
else
	error('UFp outside [0,1]')
end

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

%-Image dimensions
%-----------------------------------------------------------------------
xdim    = DIM(1); ydim = DIM(2); zdim = DIM(3);
YNaNrep = spm_type(VY(1,1).dim(4),'nanrep');

%-Intialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM    = struct(		'fname',	'mask.img',...
			'dim',		[DIM',spm_type('uint8')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM    = spm_create_image(VM);


%-Intialise beta image files
%-----------------------------------------------------------------------
Vbeta(1:nBeta) = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('float')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:nBeta
	Vbeta(i).fname   = sprintf('beta_%04d.img',i);
	Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.Xnames{i});
	spm_unlink(Vbeta(i).fname)
	Vbeta(i)         = spm_create_image(Vbeta(i));
end


%-Intialise residual sum of squares image file
%-----------------------------------------------------------------------
VResMS = struct(	'fname',	'ResMS.img',...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:Residual sum-of-squares');
VResMS = spm_create_image(VResMS);


%-Intialise residual images
%-----------------------------------------------------------------------
VResI(1:nSres) = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:Residual image'));

for i = 1:nSres
	VResI(i).fname   = sprintf('ResI_%04d.img', i);
	VResI(i).descrip = sprintf('spm_spm:ResI (%04d)', i);
	spm_unlink(VResI(i).fname);
	VResI(i)         = spm_create_image(VResI(i));
end

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...initialised')        %-#


%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Find a suitable block size for the main loop
%-----------------------------------------------------------------------
if ydim < 2, error('ydim < 2'), end		    %-need at least 2 lines
blksz  = ceil(maxMem/8/nScan/nVar);		    %-block size
nbch   = ceil(xdim*ydim/blksz);			    %-# blocks


%-Intialise other variables used through the loop 
%=======================================================================
xords  = [1:xdim]'*ones(1,ydim); xords = xords(:)'; % plane X coordinates
yords  = ones(xdim,1)*[1:ydim];  yords = yords(:)'; % plane Y coordinates
S      = 0;                                         % Volume (voxels)
Cy     = 0;					    % <Y*Y'> spatially whitened
CY     = 0;					    % <Y*Y'> for ReML
EY     = 0;					    % <Y>    for ReML
i_res  = round(linspace(1,nScan,nSres))';	    % Indices for residual

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ    = zeros(3,xdim*ydim*zdim);

%-Cycle over bunches blocks within planes to avoid memory problems
%=======================================================================
spm_progress_bar('Init',100,'model estimation','');

for z = 1:zdim				%-loop over planes (2D or 3D data)

    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords   = z*ones(xdim*ydim,1)';	%-plane Z coordinates
    Q       = {};			%-voxels analyzed (indices)
    CrBl    = [];			%-parameter estimates
    CrResI  = [];			%-normalized residuals
    CrResSS = [];			%-residual sum of squares

    for bch = 1:nbch			%-loop over blocks

	%-# Print progress information in command window
	%---------------------------------------------------------------
	fprintf('\r%-40s: %30s',sprintf('Plane %3d/%-3d, block %3d/%-3d',...
		z,zdim,bch,nbch),' ')                                %-#

	%-construct list of voxels in this block
	%---------------------------------------------------------------
	I     = [1:blksz] + (bch - 1)*blksz;		%-voxel indices
	I     = I(I <= xdim*ydim);
	xyz   = [xords(I); yords(I); zords(I)];		%-voxel coords in bch
	nVox  = size(xyz,2);				%-number of voxels

	%-Get data & construct analysis mask
	%===============================================================
	fprintf('%s%30s',sprintf('\b')*ones(1,30),'...read & mask data')%-#
	CrLm  = logical(ones(1,nVox));			%-current block mask


	%-Compute explicit mask
	% (note that these may not have same orientations)
	%---------------------------------------------------------------
	for i = 1:length(xM.VM)
		tM   = inv(xM.VM(i).mat)*M;		%-Reorientation matrix
		tmp  = tM * [xyz;ones(1,nVox)];		%-Coords in mask image

		%-Load mask image within current mask & update mask
		%-------------------------------------------------------
		CrLm(CrLm) = spm_sample_vol(xM.VM(i),...
				tmp(1,CrLm),tmp(2,CrLm),tmp(3,CrLm),0) > 0;
	end
	
	%-Get the data in mask, compute threshold & implicit masks
	%---------------------------------------------------------------
	Y     = zeros(nScan,nVox);
	for i = 1:nScan
		if ~any(CrLm), break, end		%-Break if empty mask
		Y(i,CrLm)  = spm_sample_vol(VY(i,1),... %-Load data in mask
				xyz(1,CrLm),xyz(2,CrLm),xyz(3,CrLm),0);
		CrLm(CrLm) = Y(i,CrLm) > xM.TH(i,1);	%-Threshold (& NaN) mask
		if xM.I & ~YNaNrep & xM.TH(i,1) < 0	%-Use implicit 0 mask
			CrLm(CrLm) = abs(Y(i,CrLm)) > eps;
		end
	end

	%-Mask out voxels where data is constant
	%---------------------------------------------------------------
	CrLm(CrLm) = any(diff(Y(:,CrLm),1));
	Y          = Y(:,CrLm);			%-Data matrix within mask
	CrS        = sum(CrLm);			%-# current voxels


	%===============================================================
	%-Proceed with General Linear Model (if there are voxels)
	%===============================================================
	if CrS

		%-Remove filter confounds
		%-------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...temporal filtering')     %-#

		KY         = spm_filter(xX.K,Y);

		%-General linear model: least squares estimation
		% (Using pinv to allow for non-unique designs
		% (Including temporal convolution of design matrix with K
		%------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...parameter estimation')   %-#

		beta       = xX.pKX * KY;		%-Parameter estimates
		res        = spm_sp('r',xX.xKXs,KY);	%-Residuals
		ResSS      = sum(res.^2);		%-Residual SSQ
		clear KY				%-Clear to save memory

		%-sample covariance and mean of Y (all searched voxels)
		%-------------------------------------------------------
		CY = CY + Y*Y';
		EY = EY + sum(Y,2);


		%-If UFp > 0, save raw data in 8bit squashed *.mad format
		%-------------------------------------------------------
		if UFp > 0
		fprintf('%s%30s',sprintf('\b')*ones(1,30),'saving data') %-#

			% F-threshold
			%-----------------------------------------------
			if UFp == 1			%-Save all data
			    tmp = [1:CrS];
			else
			    tmp = (sum((h*beta).^2,1)/trMV) > UF*ResSS/trRV;
			    tmp = find(tmp);
			end

			% append data and save indices to coords
			%-----------------------------------------------
			spm_append(Y(:,tmp),'Y.mad',2);
			Yidx = [Yidx, (S + tmp)];

			%-Assemble spatially whitened <Y*Y'> for ReML
			%-----------------------------------------------
			if length(tmp)
				q  = size(tmp,2);
				q  = spdiags(sqrt(trRV./ResSS(tmp)'),0,q,q);
			    	Y  = Y(:,tmp)*q;
				Cy = Cy + Y*Y';
			end

		end
		clear Y				%-Clear to save memory

		%-Save betas etc. for current plane in memory as we go along
		%-------------------------------------------------------
		CrBl 	           = [CrBl,    beta];
		CrResI             = [CrResI,  res(i_res,:)];
		CrResSS            = [CrResSS, ResSS];

	end % (CrS)

	%-Append new inmask voxel locations and volumes
	%---------------------------------------------------------------
	XYZ(:,S + [1:CrS]) = xyz(:,CrLm);	%-InMask XYZ voxel coords
	Q                  = [Q I(CrLm)];	%-InMask XYZ voxel indices
	S                  = S + CrS;		%-Volume analysed (voxels)

    end   % (bch)


    %-Plane complete, write out plane data to image files
    %===================================================================
    fprintf('%s%30s',sprintf('\b')*ones(1,30),'...saving plane')     %-#

    %-Write Mask image
    %-------------------------------------------------------------------
    tmp   = sparse(xdim,ydim);
    if length(Q)
	tmp(Q) = 1;
    end
    VM    = spm_write_plane(VM, tmp, z);

    %-Write beta images
    %-------------------------------------------------------------------
    tmp   = NaN*ones(xdim,ydim);
    for i = 1:nBeta
        if length(Q)
		tmp(Q) = CrBl(i,:);
	end
	Vbeta(i) = spm_write_plane(Vbeta(i),tmp,z);
    end

    %-Write residual images
    %-------------------------------------------------------------------
    for i = 1:nSres
	if length(Q)
		tmp(Q) = CrResI(i,:);
	end
	VResI(i) = spm_write_plane(VResI(i), tmp, z);
    end

    %-Write ResSS into ResMS (variance) image
    % (Scaling of ResSS to ResMS by tr(RV) is accomplished by adjusting the
    % (scalefactor, after the error correlations (Vi) have been estimated.
    %-------------------------------------------------------------------
    if length(Q)
	tmp(Q) = CrResSS;
    end
    VResMS     = spm_write_plane(VResMS,tmp,z);		

    %-Report progress
    %-------------------------------------------------------------------
    fprintf('%s%30s',sprintf('\b')*ones(1,30),'...done')             %-#
    spm_progress_bar('Set',100*(bch + nbch*(z - 1))/(nbch*zdim));


end % (for z = 1:zdim)
fprintf('\n')                                                        %-#
spm_progress_bar('Clear')

%=======================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%=======================================================================
if S == 0, warning('No inmask voxels - empty analysis!'), end


%-Non-sphericity: Vi
%=======================================================================
if iscell(xX.xVi.Vi)

	%-REML estimate of residual correlations through hyperparameters (h)
	%---------------------------------------------------------------
	fprintf('%-40s: %30s\n','Non-sphericity','...REML estimation') %-#
	Cy            = Cy/length(Yidx);

	% ReML for separable designs and covariance components
	%---------------------------------------------------------------
	if iscell(xX.K)
		m     = length(xX.xVi.Vi);
		h     = zeros(m,1);
		Vi    = sparse(nScan,nScan); 
		for i = 1:length(xX.K)

			% extract blocks from bases
			%-----------------------------------------------
			q     = xX.K{i}.row;
			p     = [];
			Qp    = {};
			for j = 1:m
				if any(xX.xVi.Vi{j}(q,q))
					Qp{end + 1} = xX.xVi.Vi{j}(q,q);
					p           = [p j];
				end
			end

			% design space for ReML (with confounds in filter)	
			%-----------------------------------------------
			Xp    = xX.X(q,:);
			if isfield(xX.K{i}.KH)
				Xp = [Xp xX.K{i}.KH];
			end

			% ReML
			%-----------------------------------------------
			fprintf('%-30s- %i\n','  ReML Block',i);
			[Vip,hp] = spm_reml(Cy(q,q),Xp,Qp);
			Vi(q,q)  = Vi(q,q) + Vip;
			h(p)     = hp;
		end
	else
		[Vi,h] = spm_reml(Cy,xX.X,xX.xVi.Vi);
	end

	% normalize non-sphericity and save hyperparameters
	%---------------------------------------------------------------
	Vi           = Vi*nScan/trace(Vi);
	xX.xVi.Param = h;

end


%-[Re]-enter Vi & derived values into design structure xX
%-----------------------------------------------------------------------
KVi      = spm_filter(xX.K, Vi);
xX.V     = spm_filter(xX.K,KVi'); 			%-Non-sphericity V
xX.pKXV  = xX.pKX*xX.V;					%-for contrast variance 
xX.Bcov  = xX.pKXV*xX.pKX';				%-for Cov(Beta)
[xX.trRV xX.trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);	%-Variance expectations
xX.erdf  = xX.trRV^2/xX.trRVRV;				%-Effective d.f.


%-average sample covariance and mean of Y (over voxels)
%-----------------------------------------------------------------------
CY       = CY/S;
EY       = EY/S;
CY       = CY - EY*EY';


%-Compute scaled design matrix for display purposes
%-----------------------------------------------------------------------
xX.nKX   = spm_DesMtx('sca',xX.xKXs.X,xX.Xnames);


%-Set VResMS scalefactor as 1/trRV (raw voxel data is ResSS)
%-----------------------------------------------------------------------
VResMS.pinfo(1) = 1/xX.trRV;


%-"close" written image files, updating scalefactor information
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...closing image files')  %-#
VM                      = spm_create_image(VM);
for i=1:nBeta, Vbeta(i) = spm_create_image(Vbeta(i)); end
for i=1:nSres, VResI(i) = spm_create_image(VResI(i)); end
VResMS                  = spm_create_image(VResMS);


%-Smoothness estimates of component fields and RESEL counts for volume
%=======================================================================
FWHM   = spm_est_smoothness(VResI,VM);
R      = spm_resels_vol(VM,FWHM)';

%-Delete the residuals images
%=======================================================================
for  i = 1:nSres,
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.img']);
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.hdr']);
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.mat']);
end;
clear VResI


%-Retain relative filenames for robustness
%-----------------------------------------------------------------------
fprintf('%s%30s',sprintf('\b')*ones(1,30),'...tidying file handles') %-#
VM     = VM.fname;
Vbeta  = {Vbeta.fname}';
VResMS = VResMS.fname;

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#

%-Save coordinates of within mask voxels (for Y.mad pointlist use)
%-----------------------------------------------------------------------
if UFp > 0, save Yidx Yidx, end
XYZ    = XYZ(:,1:S);

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
SPMvars = {	'SPMid','VY','xX','xM',...		%-Design parameters
		'M','DIM',...				%-Image space parameters
		'VM','Vbeta','VResMS',...		%-Filenames
		'XYZ',...				%-InMask XYZ coords
		'F_iX0','UFp','UF',...			%-F-thresholding data
		'Cy','CY',...				%-voxel-wide covariances
		'S','R','FWHM'};			%-Smoothness data

if nargin > 4

	%-Extra arguments were passed for saving in the SPM.mat file
	%---------------------------------------------------------------
	for i=5:nargin
		SPMvars = [SPMvars,{inputname(i)}];
		eval([inputname(i),' = varargin{i-4};'])
	end
end
save('SPM',SPMvars{:})


%-Save contrast structure
%-----------------------------------------------------------------------
save('xCon.mat','xCon')


fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...use the results section for assessment\n\n')             %-#
