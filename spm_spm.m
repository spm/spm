function spm_spm(VY,xX,xM,F_iX0,varargin)
% Estimation of the General Linear Model
% FORMAT spm_spm(VY,xX,xM,F_iX0,<extra parameters for SPM.mat>)
%
% VY    - Vector of structures of mapped image volumes
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
%
% xM    - Structure containing masking information, or a simple column vector
%         of thresholds corresponding to the images in VY.
%       - If a structure, the required fields are:
%         xM.TH - nScan x 1 vector of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1=>implicit zero/NaN mask)
%         xM.VM - struct array of mapped explicit mask image volumes
% 		- (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or data
%                 type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%
% F_iX0 - Indicies of design matrix columns to form the reduced design
%         matrix. The ensuing F-test (implemented through F-contrasts)
%         is labelled "effects of interest", and is to filter the saved
%         data pointlist (Y.mad).
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
% are assummed to be independent and identically Normally distributed
% with zero mean and unspecified variance. The parameters are estimated
% by ordinary least squares, using matrix methods, and the residual
% variance estimated.
%
% Note that only a single variance component is considered. This is
% *always* a fixed effects model. Random effects models can be effected
% using summary statistics, summarising the data such that the fixed
% effects in the model on the summary data are those effects of
% interest (i.e. the population effects). In this case the residual
% variance for the model on the summary data contains the appropriate
% variance components in the appropriate ratio. See spm_RandFX.man for
% further details.
%
% Under the additional assumption that the standardised residual images
% are strictly stationary standard Gaussian random fields, results from
% Random field theory can be applied to estimate the significance
% statistic images (SPM's) accounting for the multiplicity of testing
% at all voxels in the volume of interest simultaneously. The required
% parameters are the volume, and Lambda, the variance covariance matrix
% of partial derivatives of the standardised error fields. The
% algorithm estimates the variances of the partial derivatives in the
% axis directions (the diagonal of Lambda). The covariances (off
% diagonal elements of Lambda) are assumed zero.
% 
%                           ----------------
%
% For temporally correlated (fMRI) data, the algorithm implements the
% technique of Worsley & Friston (1995). In this approach the model is
% assummed to fit such that the residuals have (slight) short-term
% (intrinsic) autocorrelation given by Vi. (I.e. Residuals e =
% sqrt(Vi)*e', where e' is an (unobserved) white noise time series).
% The data and model are then temporally filtered prior to model
% fitting by a linear filter given by matrix K, giving model K*Y =
% K*X*beta + K*e. K=inv(sqrt(Vi)) corresponds to pre-whitening (leaving
% residuals e'), K=I to no temporal smoothing (appropriate for
% independent data when Vi=I), and K=S (a temporal filter matrix - see
% spm_make_filter) for band-pass filtering.
% 
% The autocorrelation in the filtered time series is then K*Vi*K'
% Standard results for serially correlated regression can then be used
% to produce variance estimates and effective degrees of freedom
% corrected for the temporal auto-correlation. Note that these
% corrections give non-integer degrees of freedom.
% 
% Vi can be passed as a parameter. If Vi is not specified then it is
% taken as I. Often Vi is unknown (or spatially variant, or just
% difficult to estimate). Provided the intrinsic autocorrelation Vi is
% small relative to that imposed by temporal smoothing (i.e. the
% temporal smoothing swamps the intrinsic autocorrelation), the
% variance-covariance structure of the residuals is approximately that
% imposed by the smoothing, i.e K*Ve*K' ~ K*Vi*K' (Ve = estimated Vi).
% if xX.xVi.row{n} is specified Vi is estimated using an AR(1) model
% for each session (n). 
% 
% The volume analsed is the intersection of the threshold masks,
% explicit masks and implicit masks. See spm_spm_ui for further details
% on masking options.
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
%     xX.V      - V matrix (xX.K*xX.xVi.Vi*xX.K')
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
%     Vbeta     - array of handle structures for beta images
%     VResMS    - handle structure of ResMS image
%     XYZ       - 3xS (S below) vector of in-mask voxel coordinates
%     F_iX0     - as input
%     UFp       - critical p-value for F-thresholding
%     UF        - F-threshold value
%     S         - Lebesgue measure or volume (in voxels)
%     R         - vector of resel counts     (in resels)
%     Lambda    - Variance-covariance matrix of partial derivatives of
%                 standardised residuals.
%                 (Covariances (off-diagonal) assummed zero)
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
%     xCon      - Contrast structure (created by spm_FcUtil.m)
%     xCon.name - Name of contrast
%     xCon.STAT - 'F' - for the F-contrast ('T' for normal contrasts)
%     xCon.c    - (F) Contrast weights
%     xCon.X0   - Reduced design matrix (spans design space under Ho)
%     xCon.iX0  - Indicies of design matrix columns to form the reduced
%                 design matrix. (Input argument F_iX0 in this case.)
%               - (Is 0 if X0 was specified, [] if c was specified.)
%     xCon.X1o  - Remaining design space (orthogonal to X0).
%     xCon.eidf - Effective interest degrees of freedom (numerator df)
%     xCon.Vcon - ...for handle of contrast/ESS image (empty at this stage)
%     xCon.Vspm - ...for handle of SPM image (empty at this stage)
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
% ``Analysis of fMRI Time-Series Revisited --- Again''
%       Neuroimage 2:173-181
%
%-----------------------------------------------------------------------
%
% For those interested, the analysis proceeds a "plank" at a time,
% where a "plank" is a bunch of lines (parallel to the x-axis) within a
% plane. The plank width is determined at run-time as the maximum
% possible (at most the plane width) plank width satisfying the memory
% constraint parameterised in the code in variable maxMem. maxMem can
% be tweaked for your system by editing the code.
%
%_______________________________________________________________________
% %W% Andrew Holmes, Jean-Baptiste Poline, Karl Friston %E%
SCCSid   = '%I%';

%-Say hello
%-----------------------------------------------------------------------
SPMid    = spm('FnBanner',mfilename,SCCSid);
Finter   = spm('FigName','Stats: estimation...'); spm('Pointer','Watch')

%-Parameters
%-----------------------------------------------------------------------
Def_UFp  = 0.001;	%-Default F-threshold for Y.mad pointlist filtering
maxMem   = 2^20;	%-Max data chunking size, in bytes
maxRes4S = 64;		%-Maximum #res images for smoothness estimation

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
				'Form',	'none'); 
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
        'Current directory already contains SPM results!',...
        '(i.e. a SPM.mat file)',...
        ['        (pwd = ',pwd,')'],...
        ' ',...
        'Existing results are being overwritten!'},...
        mfilename,sqrt(-1));
    warning(sprintf('Overwriting existing results\n\t (pwd = %s) ',pwd))
    drawnow
end

files = {	'SPM.mat','Yidx.mat','Y.mad','xCon.mat',...
		'mask.???','ResMS.???',...
		'beta_????.???','con_????.???',...
		'ess_????.???','spm?_????.???'};
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

%-Check Y images have same dimensions, orientation & voxel size
%-----------------------------------------------------------------------
if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0]) % NB: Bombs for single image
	error('images do not all have the same dimensions'), end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end



%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Initialise design space
%=======================================================================
fprintf('%-40s: %30s','Initialising design space','...checking')     %-#

%-Construct design parmameters, and store in design structure xX
% Take care to apply temporal convolution - KX stored as xX.xKX.X
%-----------------------------------------------------------------------
fprintf('%s%30s',sprintf('\b')*ones(1,30),'...computing')            %-#

[nScan nbeta] = size(xX.X);			%-#scans & #parameters

KVi           = spm_filter('apply',xX.K, xX.xVi.Vi);
xX.V          = spm_filter('apply',xX.K,KVi');  %-V matrix
KX            = spm_filter('apply',xX.K, xX.X);
xX.xKXs       = spm_sp('Set',KX);		%-Design space structure
xX.pKX        = spm_sp('x-',xX.xKXs);		%-Pseudoinverse of KX
[xX.trRV xX.trRVRV] ...				%-Variance expectations
              = spm_SpUtil('trRV',xX.xKXs,xX.V);%-(trRV & trRV2 in spm_AnCova)
xX.erdf       = xX.trRV^2/xX.trRVRV;		%-Effective residual d.f.
xX.nKX        = spm_DesMtx('sca',xX.xKXs.X,xX.Xnames); %-Scale for display

%-Check estimability
%-----------------------------------------------------------------------
if xX.erdf<0
    error(sprintf('This design is completely unestimable! (df=%-.2g)',xX.erdf))
elseif xX.erdf==0
    error('This design has no residuals! (df=0)')
elseif xX.erdf<4
    warning(sprintf('Very low degrees of freedom (df=%-.2g)',xX.erdf))
end


%-Default F-contrast (in contrast structure) & Y.mad pointlist filtering
%=======================================================================
fprintf('%s%30s',sprintf('\b')*ones(1,30),'...F-contrast')           %-#
UFp   = spm('GetGlobal','UFp'); if isempty(UFp), UFp = Def_UFp; end

if isempty(F_iX0)
	str   = 'all effects';
else
	str   = 'effects of interest';
end
xCon          = spm_FcUtil('Set',str,'F','iX0',F_iX0,xX.xKXs);
[trMV,trMVMV] = spm_SpUtil('trMV',xCon.X1o,xX.V);
xCon.eidf     = trMV^2/trMVMV;
h             = spm_FcUtil('Hsqr',xCon,xX.xKXs);

%-Compute UF, the F-threshold
%-----------------------------------------------------------------------
if UFp > 0 & UFp < 1				%-F-filter for Y.mad file
	UF   = spm_invFcdf(1 - UFp,[xCon.eidf,xX.erdf]);
	Yidx = [];
elseif UFp == 1					%-No filtering - save all data
	UF   = -Inf;
	Yidx = [];
elseif UFp == 0					%-Write no Y.mad data
	UF   = Inf;
else
	error('UFp outside [0,1]')
end

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

%-Image dimensions
%-----------------------------------------------------------------------
xdim    = VY(1).dim(1); ydim = VY(1).dim(2); zdim = VY(1).dim(3);
YNaNrep = spm_type(VY(1).dim(4),'nanrep');

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ = [];

%-Intialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM = struct(		'fname',	'mask',...
			'dim',		[VY(1).dim(1:3),spm_type('uint8')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM = spm_create_image(VM);


%-Intialise beta image files
%-----------------------------------------------------------------------
Vbeta(1:nbeta) = deal(struct(...
			'fname',	[],...
			'dim',		[VY(1).dim(1:3) spm_type('float')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i=1:nbeta
	Vbeta(i).fname   = sprintf('beta_%04d.img',i);
	Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.Xnames{i});
	spm_unlink(Vbeta(i).fname)
	Vbeta(i)         = spm_create_image(Vbeta(i));
end


%-Intialise residual sum of squares image file
%-----------------------------------------------------------------------
VResMS = struct(	'fname',	'ResMS.img',...
			'dim',		[VY(1).dim(1:3) spm_type('double')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:Residual sum-of-squares');
VResMS = spm_create_image(VResMS);

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...initialised')        %-#


%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Find a suitable block size for the main loop, which proceeds a bunch
% of lines at a time (minimum = one line; maximum = one plane)
% (maxMem is the maximum amount of data that will be processed at a time)
%-----------------------------------------------------------------------
blksz	= maxMem/8/nScan;			%-block size (in bytes)
if ydim < 2, error('ydim < 2'), end		%-need at least 2 lines
nl 	= max(min(round(blksz/xdim),ydim),1); 	%-max # lines / block
clines	= 1:nl:ydim;				%-bunch start line #'s
blines  = diff([clines ydim+1]);		%- #lines per bunch
nbch    = length(clines);			%- #bunches


%-Intialise other variables used through the loop 
%=======================================================================
BePm	= zeros(1,xdim*ydim);			%-below plane (mask)
vox_ind = [];					%-voxels indices in plane
tmp_msk	= [];					%-temporary mask 

BeVox(nbch) = struct('res',[],'ofs',[],'ind',[]);   %-voxels below
CrVox       = struct('res',[],'ofs',[],'ind',[]);   %-current voxels
						    % res : residuals
   						    % ofs : mask offset
						    % ind : indices

xords  = [1:xdim]'*ones(1,ydim); xords = xords(:)'; %-plane X coordinates
yords  = ones(xdim,1)*[1:ydim];  yords = yords(:)'; %-plane Y coordinates


%-Smoothness estimation variables
%-----------------------------------------------------------------------
S      = 0;                                     %-Volume analyzed (in voxels)
sx_res = 0; sy_res = 0; sz_res = 0;		%-sum((dr./d{x,y,z}).^2)
nx     = 0; ny     = 0; nz     = 0;		%-# {x,y,z} partial derivs

%-Indexes of residual images to sample for smoothness estimation
%-----------------------------------------------------------------------
i_res  = round(linspace(1,nScan,min(nScan,maxRes4S)))';

%-parameter for estimation of intrinsic correlations AR(1) model
%-----------------------------------------------------------------------
A      = 0;					%-regression coeficient	

%-Cycle over bunches of lines within planes (planks) to avoid memory problems
%=======================================================================
spm_progress_bar('Init',100,'model estimation','');

for z = 1:zdim				%-loop over planes (2D or 3D data)

    zords   = z*ones(xdim*ydim,1)';	%-plane Z coordinates
    CrBl    = [];			%-current plane betas
    CrResMS = [];			%-current plane ResMS
    
    for bch = 1:nbch			%-loop over bunches of lines (planks)

	%-# Print progress information in command window
	%---------------------------------------------------------------
	fprintf('\r%-40s: %30s',sprintf('Plane %3d/%-3d, plank %3d/%-3d',...
		z,zdim,bch,nbch),' ')                                %-#

	cl    = clines(bch); 	 	%-line index of first line of bunch
	bl    = blines(bch);  		%-number of lines for this bunch

	%-construct list of voxels in this bunch of lines
	%---------------------------------------------------------------
	I     = ((cl-1)*xdim+1):((cl+bl-1)*xdim);	%-lines cl:cl+bl-1
	xyz   = [xords(I); yords(I); zords(I)];		%-voxel coords in bch


	%-Get data & construct analysis mask for this bunch of lines
	%===============================================================
	fprintf('%s%30s',sprintf('\b')*ones(1,30),'...read & mask data')%-#
	CrLm  = logical(ones(1,xdim*bl));		%-current lines mask

	%-Compute explicit mask
	% (note that these may not have same orientations)
	%---------------------------------------------------------------
	for i = 1:length(xM.VM)
		tM   = inv(xM.VM(i).mat)*VY(1).mat;	%-Reorientation matrix
		tmp  = tM * [xyz;ones(1,size(xyz,2))];	%-Coords in mask image

		%-Load mask image within current mask & update mask
		%-------------------------------------------------------
		CrLm(CrLm) = spm_sample_vol(xM.VM(i),...
				tmp(1,CrLm),tmp(1,CrLm),tmp(3,CrLm),0) > 0;
	end
	
	%-Get the data in mask, compute threshold & implicit masks
	%---------------------------------------------------------------
	Y     = zeros(nScan,xdim*bl);
	for j = 1:nScan
		if ~any(CrLm), break, end		%-Break if empty mask
		Y(j,CrLm)  = spm_sample_vol(VY(j),...	%-Load data in mask
				xyz(1,CrLm),xyz(2,CrLm),xyz(3,CrLm),0);
		CrLm(CrLm) = Y(j,CrLm) > xM.TH(j);	%-Threshold (& NaN) mask
		if xM.I & ~YNaNrep & xM.TH(j)<0		%-Use implicit 0 mask
			CrLm(CrLm) = abs(Y(j,CrLm))>eps;
		end
	end

	%-Apply mask
	%---------------------------------------------------------------
	Y         = Y(:,CrLm);			%-Data matrix within mask
	CrS       = sum(CrLm);			%-#current voxels
	CrVox.ofs = I(1) - 1;			%-Current voxels line offset
	CrVox.ind = find(CrLm);			%-Voxel indicies (within bunch)



	%-Proceed with General Linear Model & smoothness estimation
	%===============================================================
	if any(CrLm)

		%-Estimate intrinsic correlation structure AR(1) model
		%-------------------------------------------------------
		switch xX.xVi.Form

		    case 'AR(1)'
		    %---------------------------------------------------
		    fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...AR(1) estimation')	     %-#

		    for i = 1:length(xX.xVi.row)
			y   = spm_detrend(Y(xX.xVi.row{i},:));
			q   = 1:(size(y,1) - 1);
			a   = sum(y(q,:).*y(q + 1,:))./sum(y(q,:).*y(q,:));
			A   = A + [1; -mean(a)];
		    end
		end

		%-Temporal smoothing
		%-------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...temporal smoothing')     %-#

		KY    = spm_filter('apply',xX.K, Y);

		%-General linear model: least squares estimation
		% (Using pinv to allow for non-unique designs            )
		% (Including temporal convolution of design matrix with K)
		%-------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...parameter estimation')   %-#
		beta  = xX.pKX * KY;			%-Parameter estimates
		res   = spm_sp('r',xX.xKXs,KY);	%-Residuals
		ResMS = sum(res.^2)/xX.trRV;		%-Residual mean square
		clear KY				%-Clear to save memory


		%-If UFp>0, save raw data in 8bit squashed *.mad file format
		%-------------------------------------------------------
		if UFp > 0
			fprintf('%s%30s',sprintf('\b')*ones(1,30),...
						'...saving data')    %-#
			if UFp == 1			%-Save all data
				spm_append(Y,'Y.mad',2)	%-append data
				Yidx = [Yidx,S+[1:CrS]];%-save indexes to coords
			else				%-F-threshold
				tmp  = (sum((h*beta).^2)/trMV) > ResMS*UF;
				spm_append(Y(:,tmp),'Y.mad',2);
				Yidx = [Yidx, (S + find(tmp))];
			end
		end
		clear Y					%-Clear to save memory


		%-Save betas for current plane in memory as we go along
		% (if there is not enough memory, could save directly as float)
		% (analyze images, or *.mad and retreived when plane complete )
		%-------------------------------------------------------
		CrBl 	= [CrBl,beta];
		CrResMS = [CrResMS,ResMS];

		% Smoothness estimation - Normalize subsampled residuals
		%-------------------------------------------------------
		res     = res(i_res,:);	
		tResSS  = sqrt(sum(res.^2));		%-ResSS of subsample
		for   j = 1:size(res,1)
			res(j,:) = res(j,:)./tResSS;
		end
		CrVox.res = res;


	   	%-Smoothness estimation: compute spatial derivatives...
		%=======================================================
		%-The code avoids looping over the voxels (using arrays)
		% and optimizes memory usage by working on the masks.
		%-It works on all the voxels contained in the mask and i_res.
		%-It constructs a list of indices that correspond
		% to the voxels in the original mask AND the 
		% displaced mask. It then finds the location of these
		% in the original list of indices. If necessary, it 
		% constructs the voxel list corresponding to the displaced 
		% mask (eg in y-dim when not at the begining of the plane).

		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...smoothness estimation') %-#


		%-Construct utility vector to map from mask (CrLm) to voxel
		% index within the mask. (For voxels with CrLm true, j_jk
		% is the index of the voxel within the mask)
		j_jk = cumsum(CrLm);			% (valid for x y z)


		%- z dim
		%-------------------------------------------------------
		if z>1			%-There's a plane below
		    %-Indices of inmask voxels with inmask z-1 neighbours
		    jk = find(CrLm & BePm(I));
		    if ~isempty(jk)
			%-Compute indices of inmask z-1 voxels
			k_jk   = cumsum(BePm(I)); 
			sz_res = sz_res + sum(sum( ...
					(CrVox.res(:,j_jk(jk)) ...
					- BeVox(bch).res(:,k_jk(jk))).^2 ));
			nz     = nz + length(jk);
		   end % (if ~isempty(jk))
		end % (if z~=1)


		%- y dim
		%-------------------------------------------------------
		if bch==1		%-first bunch : no previous line
		    if bl > 1
			%-Indices of inmask voxels with inmask y+1 neighbours
		   	jk = find(CrLm & [CrLm(xdim+1:xdim*bl),zeros(1,xdim)]);
		   	if ~isempty(jk)
			    %-Compute indices of inmask y+1 voxels
			    k_jk   = cumsum(CrLm(xdim+1:xdim*bl))+j_jk(xdim);
			    % NB: sum(CrLm(1:xdim))==j_jk(xdim)
			    %-Compute partial derivs, sum ^2's over resids & vox
			    sy_res = sy_res + sum(sum( ...
				(CrVox.res(:,k_jk(jk)) - ...
			 	 CrVox.res(:,j_jk(jk)) ).^2  ));
		      	   ny      = ny + length(jk);
		   	end % (if ~isempty(jk))
		    end % (if bl > 1)


		else % (if bch==1)	%-a previous line exists
		    %-Make tmp_mask as mask shifted y-1, by prepending previous
		    % line to CrLm minus its last line
		    tmp_msk = [BePm( ((cl-2)*xdim+1):(cl-1)*xdim ),...
				CrLm(1:xdim*(bl-1)) ];

		    %-get residuals for y-1 shifted space (of tmp_vox)
		    tmp_vox = ...
		    [BeVox(bch-1).res(:,find(BeVox(bch-1).ind>xdim*(nl-1))),...
		       CrVox.res(:,find(CrVox.ind <= xdim*(bl-1))) ];
	
		    %-Inmask voxels with inmask (incl. prev line) y-1 neighbours
		    jk	= find(CrLm & tmp_msk);
		    if ~isempty(jk)
			%-Compute indices of inmask y-1 voxels
			k_jk	= cumsum(tmp_msk);
			sy_res 	= sy_res + sum(sum( ...
					(CrVox.res(:,j_jk(jk)) - ...
					tmp_vox(:,k_jk(jk))).^2  ));
			ny 	= ny + length(jk);
		    end % (if ~isempty(jk))

		end % (if bch==1)


		%- x dim
		%-------------------------------------------------------
		%-Shift the mask to the left, add 0 at line ends
		tmp_msk = [CrLm(2:length(CrLm)),0];
		tmp_msk(xdim:xdim:bl*xdim) = 0;
		
		%-Indices inmask of voxels with inmask x+1 neighbours
		jk = find(CrLm & tmp_msk);
		if ~isempty(jk)
		   sx_res = sx_res + sum(sum( ...
		   	(CrVox.res(:,j_jk(jk)+1) - ...
			 CrVox.res(:,j_jk(jk))).^2  ));
		   % NB: j_jk(jk)+1 = position of voxels on right of j_jk(jk)
		   nx     = nx + length(jk);
		end % (if ~isempty(jk))

	    end % (if any(CrLm))

	    %-Append new inmask voxel locations and volumes
	    %-----------------------------------------------------------
	    XYZ            = [XYZ,xyz(:,CrLm)];	%-InMask XYZ voxel coordinates
	    S              = S + CrS;		%-Volume analysed (voxels)
	    					% (equals size(XYZ,2))

	    %-Roll... (BeVox(bch) is overwritten)
	    %-----------------------------------------------------------
	    BeVox(bch).ind = CrVox.ind;		%-Voxel indexes (within bunch)
	    BeVox(bch).ofs = CrVox.ofs;		%-Bunch voxel offset (in plane)
	    BeVox(bch).res = CrVox.res;		%-Sample of residuals
	    BePm(I)        = CrLm;		%-"below plane" mask

    end % (for bch = 1:nbch)
    

    %-Plane complete, write out plane data to image files
    %===================================================================
    fprintf('%s%30s',sprintf('\b')*ones(1,30),'...saving plane')     %-#

    %-Mask image
    %-BePm (& BeVox) now contains a complete plane mask
    %-------------------------------------------------------------------
    VM    = spm_write_plane(VM, reshape(BePm,xdim,ydim), z);

    %-Construct voxel indices for BePm
    %-------------------------------------------------------------------
    Q     = find(BePm);
    tmp   = NaN*ones(xdim,ydim);


    %-Write beta images
    %-------------------------------------------------------------------
    for i = 1:nbeta
        if length(Q), tmp(Q) = CrBl(i,:); end
	Vbeta(i) = spm_write_plane(Vbeta(i),tmp,z);
    end
	    
    %-Write ResMS (variance) image
    %-------------------------------------------------------------------
    if length(Q), tmp(Q) = CrResMS; end
    VResMS = spm_write_plane(VResMS,tmp,z);		   

    %-Report progress
    %-------------------------------------------------------------------
    fprintf('%s%30s',sprintf('\b')*ones(1,30),'...done')             %-#
    spm_progress_bar('Set',100*(bch + nbch*(z-1))/(nbch*zdim));

end % (for z = 1:zdim)
fprintf('\n')                                                        %-#


%=======================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%=======================================================================
if S == 0, warning('No inmask voxels - empty analysis!'), end

%-Intrinsic autocorrelations Vi (session specific)
%-----------------------------------------------------------------------
switch xX.xVi.Form

	case 'AR(1)'
	%---------------------------------------------------
	p     = length(A) - 1;			% order AR(p)
	A     = A/A(1);
	for i = 1:length(xX.xVi.row)
		q     = xX.xVi.row{i};
		n     = length(q);
		Ki    = inv(spdiags(ones(n,1)*A',-[0:p],n,n));
		Ki    = Ki.*(Ki > 1e-6);
		Vi    = Ki*Ki';
		D     = spdiags(sqrt(1./diag(Vi)),0,n,n);
		Vi    = D*Vi*D;
		xX.xVi.Vi(q,q) = Vi;
	end
end
xX.xVi.Param = A;


%-[Re]-enter Vi in design structure xX and default contrast
%-----------------------------------------------------------------------
KVi           = spm_filter('apply',xX.K, xX.xVi.Vi);
xX.V          = spm_filter('apply',xX.K,KVi');  %-V matrix
xX.pKXV       = xX.pKX*xX.V;			%-for contrast variance weight
xX.Bcov       = xX.pKXV*xX.pKX';		%-Variance of est. param.
[trRV trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);%-Variance expectations
xX.erdf       = trRV^2/trRVRV;			%-Effective residual d.f.

[trMV trMVMV] = spm_SpUtil('trMV',xCon.X1o,xX.V);
xCon.eidf     = trMV^2/trMVMV;

%-adjust trRV in ResMS scalefactor ans design structure
%-----------------------------------------------------------------------
VResMS.pinfo(1) = xX.trRV/trRV;
xX.trRV         = trRV;
xX.trRVRV       = trRVRV;


%-Smoothness estimates of component fields
%-----------------------------------------------------------------------
if zdim == 1
	if any(~[nx,ny])
		warning(sprintf('W: nx=%d, ny=%d',nx,ny)), end
	Lambda = diag([sx_res/nx sy_res/ny Inf]*(xX.erdf-2)/(xX.erdf-1));
else
	if any(~[nx,ny,nz])
		warning(sprintf('W: nx=%d, ny=%d, nz=%d',nx,ny,nz)), end
	Lambda = diag([sx_res/nx sy_res/ny sz_res/nz]*(xX.erdf-2)/(xX.erdf-1));
end
W      = (2*diag(Lambda)').^(-1/2);
FWHM   = sqrt(8*log(2))*W;


%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#

%-"close" image files, updating scalefactor information
%-----------------------------------------------------------------------
VM                      = spm_create_image(VM);
for i=1:nbeta, Vbeta(i) = spm_create_image(Vbeta(i)); end
VResMS                  = spm_create_image(VResMS);


%-Estimate RESEL counts for volume
%-----------------------------------------------------------------------
% R    = spm_resels_vol(VM,FWHM);
R      = spm_resels(FWHM,S);		%-Assumme a sphere


%-Save coordinates of within mask voxels (for Y.mad pointlist use)
%-----------------------------------------------------------------------
if UFp > 0, save Yidx Yidx, end

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
SPMvars = {	'SPMid','VY','xX','xM',...	%-General design parameters
		'Vbeta','VResMS',...		%-Handles of beta & ResMS images
		'XYZ',...			%-InMask XYZ voxel coords
		'F_iX0','UFp','UF',...		%-F-thresholding data
		'S','R','Lambda','W','FWHM'};	%-Smoothness data

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
spm_progress_bar('Clear')
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...use the results section for assessment\n\n')             %-#
