function [SPM] = spm_spm(SPM)
% Estimation of the General Linear Model
% FORMAT [SPM] = spm_spm(SPM)
%
% required fields of SPM:
%
% xY.VY - nScan x 1 struct array of mapped image volumes
%         Images must have the same orientation, voxel size and data type
%       - Any scaling should have already been applied via the image handle
%         scalefactors.
%
% xX    - Structure containing design matrix information
%       - Required fields are:
%         xX.X      - Design matrix (raw, not temporally smoothed)
%         xX.name   - cellstr of parameter names corresponding to columns
%                     of design matrix
%       - Optional fields are:
%         xX.K      - Sparse temporal smoothing matrix (see spm_sptop)
% 		    - Or cell of session-specific structures (see spm_filter)
%                   - Design & data are smoothed using K
%                     (model is K*Y = K*X*beta + K*e)
%                   - Note that K should not smooth across block boundaries
%                   - defaults to speye(size(xX.X,1))
%
% xVi   - structure describing intrinsic temporal auto-correlation
%       - required fields are
%         SPM.xVi.Vi - intrinsic temporal auto-correlation matrix
%                    - defaults to speye(size(xX.X,1))
%                    - specifying a cell array of contraints
%                      invokes spm_reml to estimate hyperparameters
%                      assuming Vi is constant over voxels.
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
% xCon  - See Christensen for details of F-contrasts
%
%
% In addition, the global default UFp is used to set a critical
% F-threshold for selecting voxels over which the non-spherical
% correlations structure is estimated (if required)
%
%_______________________________________________________________________
%
% spm_spm is the heart of the SPM package. Given image files and a
% General Linear Model, it estimates the model parameters, residual
% variance, and smoothness of the standardised residual fields, writing
% these out to disk in the current working directory for later
% interrogation in the results section. (NB: Existing analyses in the
% current working directory are overwritten.).  This directory
% now becomes the working directory for this analysis and all saved
% images are relative to this directory.
%
% The model is expressed via the design matrix (xX.X). The basic model
% at each voxel is of the form is Y = X*beta * e, for data Y, design
% matrix X, (unknown) parameters B, and residual errors e. The errors
% are assummed to have a normal distribution ~ N(0,sigma^2*Vi) where Vi can
% be specified in xVi.  If the intrinsic correlations Vi are not
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
%
% spm_est_smoothness estimates the variances of the partial derivatives 
% in the axis directions (the diagonal of Lambda). The covariances (off
% diagonal elements of Lambda) are assumed to be zero.
% 
%                           ----------------
%
% In some instances the i.i.d. assumptions about errors do not hold.  For
% example, with serially correlated (fMRI) data or correlations among the
% levels of a factor in repeated measures designs.  This non-sphericity
% can be specified in terms of constraints (SPM.xVi.Vi) Covariance
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
% The output of spm_spm takes the form of an SPM.mat file of the analysis
% parameters, and 'float' flat-file images of the parameter and variance
% [hyperparameter] estimates. An 8bit zero-one mask image indicating the
% voxels assessed is also written out, with zero indicating voxels outside
% tha analysed volume.
%
%                           ----------------
%
% The following SPM.fields are computed by spm_spm
%
%     Vbeta     - struct array of beta image handles (relative)
%     VResMS    - file struct of ResMS image handle  (relative)
%     VM        - file struct of Mask  image handle  (relative)
%
%                           ----------------
%
%     xX.V      - V matrix (xX.K*Vi*xX.K') = correlations after K is applied
%     xX.xKXs   - space structure for K*X (xX.K*xX.X), the 'filtered'
%                 design matrix
%               - given as spm_sp('Set',xX.K*xX.X) - see spm_sp for details
%                 (so xX.xKXs.X is K*X)
%     xX.pKX    - pseudoinverse of K*X (xX.K*xX.X), computed by spm_sp
%                 (this was called PH in spm_AnCova)
%     xX.Bcov   - pinv(K*X)*V*pinv(K*X)' - variance-covariance matrix
%                 of parameter estimates
%		  (when multiplied by the voxel-specific hyperparameter ResMS)
%                 (of the parameter estimates. (ResSS/xX.trRV = ResMS)    )
%     xX.trRV   - trace of R*V, computed efficiently by spm_SpUtil
%     xX.trRVRV - trace of RVRV
%     xX.erdf   - effective residual degrees of freedom (trRV^2/trRVRV)
%     xX.nKX    - design matrix (xX.xKXs.X) scaled for display
%                 (see spm_DesMtx('sca',... for details)
%
%                           ----------------
%
%     xVol.M    - 4x4 voxel->mm transformation matrix
%     xVol.iM   - 4x4 mm->voxel transformation matrix
%     xVol.DIM  - image dimensions - column vector (in voxels)
%     xVol.XYZ  - 3 x S vector of in-mask voxel coordinates
%     xVol.S    - Lebesgue measure or volume       (in voxels)
%     xVol.R    - vector of resel counts           (in resels)
%     xVol.FWHM - Smoothness of components - FWHM, (in voxels)
%
%                           ----------------
%
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
%     xVi.V      - esimated non-spherical correlation matrix
%     xVi.h      - hyperparameters  xVi.V = xVi.h(1)*xVi.Vi{1} + ...
%     xVi.Cy     - spatially whitened <Y*Y'> (used by ReML to estimate h)
%     xVi.CY     - <(Y - <Y>)*(Y - <Y>)'>    (used by spm_spm_Bayes)
%
%                           ----------------
%
% The following images are written to file
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


%-Delete files from previous analyses
%-----------------------------------------------------------------------
if exist(fullfile('.','mask.img'),'file') == 2
	spm('alert!',{...
		'Current directory contains some SPM results files!',...
    		['        (pwd = ',pwd,')'],...
		'Existing results are being overwritten!'},...
		mfilename,1);
	warning(sprintf('Overwriting existing results\n\t (pwd = %s) ',pwd))
	drawnow
end

files = {	'mask.???','ResMS.???','RVP.???',...
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


%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Initialise
%=======================================================================
fprintf('%-40s: %30s','Initialising parameters','...computing')    %-#
xX            = SPM.xX;
[nScan nBeta] = size(xX.X);


%-Check confounds (K) and non-sphericity (xVi) fields of xX structure
%-----------------------------------------------------------------------
if ~isfield(xX,'K')
	xX.K = speye(size(xX.X,1));
end
if ~isfield(SPM,'xVi')
	xVi  = struct(	'Vi',	speye(size(xX.X,1)),...
				'form',	'none'); 
else
	xVi  = SPM.xVi;
end

%-If xM is not a structure then assumme it's a vector of thresholds
%-----------------------------------------------------------------------
xM       = SPM.xM;
if ~isstruct(xM), xM = struct(	'T',	[],...
				'TH',	xM,...
				'I',	0,...
				'VM',	{[]},...
				'xs',	struct('Masking','analysis threshold'));
end


%-Image dimensions and data
%=======================================================================
VY       = SPM.xY.VY;
M        = VY(1).mat;
DIM      = VY(1).dim(1:3)';
xdim     = DIM(1); ydim = DIM(2); zdim = DIM(3);
YNaNrep  = spm_type(VY(1).dim(4),'nanrep');

global defaults
%-Maximum number of residual images for smoothness estimation
%-----------------------------------------------------------------------
MAXRES   = defaults.stats.maxres;

%-maxMem is the maximum amount of data processed at a time (bytes)
%-----------------------------------------------------------------------
MAXMEM   = defaults.stats.maxmem;

nSres    = min(nScan,MAXRES);
blksz    = ceil(MAXMEM/8/nScan);				%-block size
nbch     = ceil(xdim*ydim/blksz);				%-# blocks


%-Parameters for estimating non-sphericity (based on first F-contrast)
%-----------------------------------------------------------------------
trMV     = spm_SpUtil('trMV',spm_FcUtil('X1o',SPM.xCon(1),xX.xKXs));
trRV     = spm_SpUtil('trRV',xX.xKXs);
h        = spm_FcUtil('Hsqr',SPM.xCon(1),xX.xKXs);

%-Compute UF, the F-threshold
%-----------------------------------------------------------------------
UFp      = spm('GetGlobal','UFp');
UF       = spm_invFcdf(1 - UFp,[trMV,trRV]);

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

%-Intialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM    = struct(		'fname',	'mask.img',...
			'dim',		[DIM',spm_type('uint8')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM    = spm_create_vol(VM);


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
	Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.name{i});
	spm_unlink(Vbeta(i).fname)
end
Vbeta = spm_create_vol(Vbeta,'noopen');


%-Intialise residual sum of squares image file
%-----------------------------------------------------------------------
VResMS = struct(	'fname',	'ResMS.img',...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:Residual sum-of-squares');
VResMS = spm_create_vol(VResMS,'noopen');


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
end
VResI = spm_create_vol(VResI);

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...initialised')        %-#


%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================


%-Intialise variables used in the loop 
%=======================================================================
xords = [1:xdim]'*ones(1,ydim); xords = xords(:)';  % plane X coordinates
yords = ones(xdim,1)*[1:ydim];  yords = yords(:)';  % plane Y coordinates
S     = 0;                                          % Volume (voxels)
s     = 0;                                          % Volume (voxels > UF)
Cy    = 0;					    % <Y*Y'> spatially whitened
CY    = 0;					    % <Y*Y'> for ReML
EY    = 0;					    % <Y>    for ReML
i_res = round(linspace(1,nScan,nSres))';	    % Indices for residual

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ   = zeros(3,xdim*ydim*zdim);

%-Cycle over bunches blocks within planes to avoid memory problems
%=======================================================================
spm_progress_bar('Init',100,'model estimation','');

for z = 1:zdim				%-loop over planes (2D or 3D data)

    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords   = z*ones(xdim*ydim,1)';	%-plane Z coordinates
    CrBl    = [];			%-parameter estimates
    CrResI  = [];			%-normalized residuals
    CrResSS = [];			%-residual sum of squares
    Q       = [];			%-in mask indices for this plane

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
	Cm    = logical(ones(1,nVox));			%-current mask


	%-Compute explicit mask
	% (note that these may not have same orientations)
	%---------------------------------------------------------------
	for i = 1:length(xM.VM)

		%-Coordinates in mask image
		%-------------------------------------------------------
		tmp    = inv(xM.VM(i).mat)*M*[xyz;ones(1,nVox)];

		%-Load mask image within current mask & update mask
		%-------------------------------------------------------
		Cm(Cm) = spm_get_data(xM.VM(i),tmp(:,Cm)) > 0;
	end
	
	%-Get the data in mask, compute threshold & implicit masks
	%---------------------------------------------------------------
	Y     = zeros(nScan,nVox);
	for i = 1:nScan

		%-Load data in mask
		%-------------------------------------------------------
		if ~any(Cm), break, end			%-Break if empty mask
		Y(i,Cm) = spm_get_data(VY(i),xyz(:,Cm));

		Cm(Cm)  = Y(i,Cm) > xM.TH(i);		%-Threshold (& NaN) mask
		if xM.I & ~YNaNrep & xM.TH(i) < 0	%-Use implicit mask
			Cm(Cm) = abs(Y(i,Cm)) > eps;
		end
	end

	%-Mask out voxels where data is constant
	%---------------------------------------------------------------
	Cm(Cm) = any(diff(Y(:,Cm),1));
	Y      = Y(:,Cm);				%-Data within mask
	CrS    = sum(Cm);				%-# current voxels


	%===============================================================
	%-Proceed with General Linear Model (if there are voxels)
	%===============================================================
	if CrS

		%-Remove filter confounds
		%-------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...temporal filtering')     %-#

		KY    = spm_filter(xX.K,Y);

		%-General linear model: least squares estimation
		% Using pinv to allow for non-unique designs
		% Including temporal convolution of design matrix with K
		%------------------------------------------------------
		fprintf('%s%30s',sprintf('\b')*ones(1,30),...
					'...parameter estimation')   %-#

		beta  = xX.pKX*KY;			%-Parameter estimates
		res   = spm_sp('r',xX.xKXs,KY);		%-Residuals
		ResSS = sum(res.^2);			%-Residual SSQ
		clear KY				%-Clear to save memory


		%-sample covariance and mean of Y (all searched voxels)
		%-------------------------------------------------------
		CY    = CY + Y*Y';
		EY    = EY + sum(Y,2);


		%-F-threshold & accumulate spatially whitened Y*Y' for ReML
		%-------------------------------------------------------
		tmp   = sum((h*beta).^2,1)/trMV > UF*ResSS/trRV;
		tmp   = find(tmp);
		if length(tmp)
			q  = size(tmp,2);
			s  = s  + q;
			q  = spdiags(sqrt(trRV./ResSS(tmp)'),0,q,q);
			Y  = Y(:,tmp)*q;
			Cy = Cy + Y*Y';
		end
		clear Y				%-Clear to save memory

		%-Save betas etc. for current plane in memory as we go along
		%-------------------------------------------------------
		CrBl 	   = [CrBl,    beta];
		CrResI     = [CrResI,  res(i_res,:)];
		CrResSS    = [CrResSS, ResSS];

	end % (CrS)

	%-Append new inmask voxel locations and volumes
	%---------------------------------------------------------------
	XYZ(:,S + [1:CrS]) = xyz(:,Cm);		%-InMask XYZ voxel coords
	Q                  = [Q I(Cm)];		%-InMask XYZ voxel indices
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
	Vbeta(i) = spm_write_plane(Vbeta(i), tmp, z);
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
if iscell(xVi.Vi)

	%-REML estimate of residual correlations through hyperparameters (h)
	%---------------------------------------------------------------
	fprintf('%-40s: %30s\n','Non-sphericity','...REML estimation') %-#
	Cy            = Cy/s;

	% ReML for separable designs and covariance components
	%---------------------------------------------------------------
	if isstruct(xX.K)
		m     = length(xVi.Vi);
		h     = zeros(m,1);
		Vi    = sparse(nScan,nScan); 
		for i = 1:length(xX.K)

			% extract blocks from bases
			%-----------------------------------------------
			q     = xX.K(i).row;
			p     = [];
			Qp    = {};
			for j = 1:m
				if any(xVi.Vi{j}(q,q))
					Qp{end + 1} = xVi.Vi{j}(q,q);
					p           = [p j];
				end
			end

			% design space for ReML (with confounds in filter)	
			%-----------------------------------------------
			Xp       = xX.X(q,:);
			if isfield(xX.K(i).KH)
				Xp = [Xp xX.K(i).KH];
			end

			% ReML
			%-----------------------------------------------
			fprintf('%-30s- %i\n','  ReML Block',i);
			[Vip,hp] = spm_reml(Cy(q,q),Xp,Qp);
			Vi(q,q)  = Vi(q,q) + Vip;
			h(p)     = hp;
		end
	else
		[Vi,h] = spm_reml(Cy,xX.X,xVi.Vi);
	end

	% normalize non-sphericity and save hyperparameters
	%---------------------------------------------------------------
	Vi    = Vi*nScan/trace(Vi);
	xVi.h = h;
	xVi.V = Vi;

else
	Vi    = xVi.Vi;
	xVi.V = Vi;
end


%-[Re]-enter K*Vi*K (xX.V) & derived values into design structure xX
%-----------------------------------------------------------------------
xX.V          = spm_filter(xX.K,spm_filter(xX.K,Vi)'); 	%-Non-sphericity V
xX.Bcov       = xX.pKX*xX.V*xX.pKX';			%-Cov(Beta)
[trRV trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);	%-Expectations
xX.trRV       = trRV;					% <R'*y'*y*R>
xX.trRVRV     = trRVRV;					% for
xX.erdf       = trRV^2/trRVRV;				%-Satterthwiate approx..


%-average sample covariance and mean of Y (over voxels)
%-----------------------------------------------------------------------
CY            = CY/S;
EY            = EY/S;
CY            = CY - EY*EY';

%-Compute scaled design matrix for display purposes
%-----------------------------------------------------------------------
xX.nKX        = spm_DesMtx('sca',xX.xKXs.X,xX.name);


%-close written image files, updating scalefactor information
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...closing image files')  %-#
VM              = spm_close_vol(VM);
Vbeta           = spm_close_vol(Vbeta);
VResI           = spm_close_vol(VResI);

%-Set VResMS scalefactor as 1/trRV (raw voxel data is ResSS)
%-----------------------------------------------------------------------
VResMS.pinfo(1) = 1/xX.trRV;
VResMS          = spm_close_vol(VResMS);
VResMS          = spm_create_vol(VResMS,'noopen');
VResMS          = spm_close_vol(VResMS);


%-Smoothness estimates of component fields and RESEL counts for volume
%=======================================================================
[FWHM,VRpv]  = spm_est_smoothness(VResI,VM);
R            = spm_resels_vol(VM,FWHM)';

%-Delete the residuals images
%=======================================================================
for  i = 1:nSres,
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.img']);
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.hdr']);
	spm_unlink([spm_str_manip(VResI(i).fname,'r') '.mat']);
end;
clear VResI

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#

%-place fields in SPM
%-----------------------------------------------------------------------
SPM.xVol.XYZ   = XYZ(:,1:S);			%-InMask XYZ coords (voxels)
SPM.xVol.M     = M;				%-voxels -> mm
SPM.xVol.iM    = inv(M);			%-mm -> voxels
SPM.xVol.DIM   = DIM;				%-image dimensions
SPM.xVol.FWHM  = FWHM;				%-Smoothness data
SPM.xVol.R     = R;				%-Resel counts
SPM.xVol.S     = S;				%-Volume (voxels)

SPM.Vbeta      = Vbeta;				%-Filehandle - Beta
SPM.VResMS     = VResMS;			%-Filehandle - Hyperparameter
SPM.VM         = VM;				%-Filehandle - Mask
SPM.VRpv       = VRpv;				%-Filehandle - Resels per voxel

SPM.xVi        = xVi;				% non-sphericity structure
SPM.xVi.Cy     = Cy;				%-spatially whitened <Y*Y'>
SPM.xVi.CY     = CY;				%-<(Y - <Y>)*(Y - <Y>)'> 

SPM.xX         = xX;				%-design structure
SPM.xM         = xM;				%-mask structure

SPM.SPMid      = SPMid;
SPM.swd        = pwd;

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
save SPM SPM

%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...use the results section for assessment\n\n')             %-#
