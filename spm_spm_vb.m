function [SPM] = spm_spm_vb(SPM)
% VB estimation of voxel-specific GLM-AR models with spatial regularisation
% FORMAT [SPM] = spm_spm_vb(SPM)
%
% This function implements a VB estimation scheme for voxel-specific GLM-AR
% models. Both regression coefficients and AR coefficients are spatially
% regularised. The algorithm is described in a series of papers:
%
% Paper VB1: W. Penny, S. Kiebel and K. Friston (2003) Variational Bayesian 
%            Inference for fMRI time series. NeuroImage 19, pp 727-741.
%
% Paper VB2: W.Penny, N. Trujillo-Bareto and K. Friston (2004). Bayesian fMRI 
%            time series analysis with spatial priors. Submitted to NeuroImage.
%
% Paper VB3: W.Penny (2004). Bayesian analysis of single-subject fMRI data: 
%            SPM implementation. Technical Report. WDIN, UCL.
%
% Paper VB4: W.Penny et al. (2004) Beyond the voxel: Comparing spatially extended
%            fMRI models. Manuscript in preparation.
%
% The space to be analysed is a 'Volume', 'Masked Volume' or 'Slices'
% For 'Slices' the numbers of the slices to be analysed are then entered
%
% ________________________________________________________________________________
%
% Required fields of SPM:
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
%
% xM    - Structure containing masking information, or a simple column vector
%         of thresholds corresponding to the images in VY.
%       - If a structure, the required fields are:
%         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1 => implicit zero/NaN mask)
%         xM.VM - struct array of mapped explicit mask image volumes
% 		- (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or data
%                 type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%       - Note that voxels with constant data (i.e. the same value across
%         scans) are also automatically masked out.
%         
% ________________________________________________________________________________
%
% spm_spm_vb adds the following fields to SPM:
%
%	SPM.VCbeta     - Handles of posterior parameter estimates (Cbeta_????)
%   SPM.VPsd       - Handles of SD of posterior parameter estimates (SDbeta_????)
%	SPM.VHp        - Handle of standard deviation of the error (SDerror)
%   SPM.VAR        - Handles of AR coefficient images (AR_????)
%
%   SPM.PPM        - Posterior Probability Map data structure
%
%                   .VB=1, tells later functions (spm_contrasts, spm_graph) 
%                    that parameters were estimated using VB
%                       
%                   .AR_P, assumed AR model order
%
%                   .priors, type of priors used (eg. 'Spatial-GMRF')
%                            see spm_vb_set_priors.m
%
%                   .compute_F, whether model evidence was computed
%
%                   .Gamma, default effect size threshold (used in spm_getSPM)
%
%                   .M_X, set to 1 if an explicit mask is used (otherwise 0)
%
%                   The following parameters are set if the space
%                   to be analysed chosen as 'Slices'
%                       .AN_slices, numbers of slices analysed
%                       .slice(z), further info about GLM-AR model at slice z
%                                  eg. slice(z).F is evidence for slice z 
%                                  (if computed)
%                       These fields are used later in spm_contrasts.m, spm_graph.m
%
% ________________________________________________________________________________                         
%
% The following images are written to file:
%
% mask.{img,hdr}                                   - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
%
% Note mask.img is only written if the selected space is 'Volume' or 
% 'Masked Volume' (ie not 'Slices')
%
% Cbeta_????.{img,hdr}  
% These are 16-bit (float) images of the parameter posteriors. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.img) are given
% value NaN.
%
% SDbeta_????.{img,hdr} 
% These are 16-bit (float) images of the standard deviation of parameter 
% posteriors. 
% The image files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.img) are given
% value NaN.
%
% SDerror.{img,hdr} 
% This is a 16-bit (float) image of the standard deviation of the error.
% Voxels outside the analysis mask (mask.img) are given value NaN.
%
% AR_????.{img,hdr} 
% This is a 16-bit (float) image of the AR coefficient.
% The image files are numbered according to the order of the corresponding
% AR coefficient.
% Voxels outside the analysis mask (mask.img) are given value NaN.
%
%                           ----------------
%
% %W% Will Penny %E%

% Let later functions know (eg. spm_contrasts) that 
% estimation was with VB
PPM.VB=1;

%-Say hello
%-----------------------------------------------------------------------
Finter   = spm('FigName','Stats: Bayesian Specification ...'); 
spm('Pointer','Arrow');

%-Get SPM.mat if necessary
%-----------------------------------------------------------------------
if nargin ==0
	swd     = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');
	load(fullfile(swd,'SPM.mat'));
	SPM.swd = swd;
end

%-Change to SPM.swd if specified
%-----------------------------------------------------------------------
try
	cd(SPM.swd);
end

%-Ensure data are assigned
%-----------------------------------------------------------------------
try
	SPM.xY.VY;
catch
	helpdlg({	'Please assign data to this design',...
			'Use fMRI under model specification'});
	spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
	return
end

%-Delete files from previous analyses
%-----------------------------------------------------------------------
if exist(fullfile('.','mask.img'),'file') == 2

	str   = {'Current directory contains SPM estimation files:',...
		 'pwd = ',pwd,...
		 'Existing results will be overwritten!'};

	abort = spm_input(str,1,'bd','stop|continue',[1,0],1);
	if abort
		spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
		return
	else
		str = sprintf('Overwriting old results\n\t (pwd = %s) ',pwd)
		warning(str)
		drawnow
	end
end

files = {	'mask.???','ResMS.???','RVP.???',...
		'beta_????.???','con_????.???','ResI_????.???',...
		'ess_????.???', 'spm?_????.???'};

for i=1:length(files)
	if any(files{i} == '*'|files{i} == '?' )
		[j,null] = spm_list_files(pwd,files{i});
		for i=1:size(j,1)
			spm_unlink(deblank(j(i,:)))
		end
	else
		spm_unlink(files{i})
	end
end

%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================



%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

VY       = SPM.xY.VY;
M        = VY(1).mat;
DIM      = VY(1).dim(1:3)';
xdim   = DIM(1); ydim = DIM(2); zdim = DIM(3);

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ   = zeros(3,xdim*ydim*zdim);

%-Intialise conditional estimate image files
%-----------------------------------------------------------------------
xX             = SPM.xX;
[nScan nBeta]  = size(xX.X);
Vbeta(1:nBeta) = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('float')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:nBeta
	Vbeta(i).fname   = sprintf('Cbeta_%04d.img',i);
	Vbeta(i).descrip = sprintf('Posterior mean of beta (%04d) - %s',i,xX.name{i});
	spm_unlink(Vbeta(i).fname)
end
Vbeta = spm_create_vol(Vbeta,'noopen');

%-Intialise hyperparameter (AR 1..p and noise variance) image files
%-----------------------------------------------------------------------

% Set number of AR coefficients
str=' AR model order';
PPM.AR_P      = spm_input(str,1,'e',0,[1 1]);

% Specify type of prior
spm_input('Overall Priors ...',1,'d');
Ctype = {
    'Spatial - GMRF',...
    'Spatial - LORETA',...
        'Voxel - Shrinkage',...
        'Voxel - Uninformative'};
str   = 'Select prior';
Sel   = spm_input(str,2,'m',Ctype);
PPM.priors.WA = Ctype{Sel};

% Override option for prior on A - none
PPM.priors.overA='none';

% Compute evidence ?
str = 'Compute evidence ?';
PPM.compute_F = spm_input(str,1,'b',{'yes','no'},[1 0]);

% Set number of VB iterations
%str='Number of VB iterations';
%maxVBits      = spm_input(str,1,'e',0,[1 1]);

% Initialise Error SD image
VHp        = deal(struct(...
    'fname',	[],...
    'dim',		[DIM',spm_type('double')],...
    'mat',		M,...
    'pinfo',	[1 0 0]',...
    'descrip',	''));

VHp.fname   = sprintf('SDerror.img');
VHp.descrip = sprintf('Error SD');
spm_unlink(VHp.fname)
VHp   = spm_create_vol(VHp,'noopen');

% Initialise Posterior SD images
nPsd=size(SPM.xX.X,2);
VPsd(1:nPsd)        = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:nPsd
	VPsd(i).fname   = sprintf('SDbeta_%04d.img',i);
	VPsd(i).descrip = sprintf('Posterior SD of beta (%04d)',i);
	spm_unlink(VPsd(i).fname)
end
VPsd   = spm_create_vol(VPsd,'noopen');

% Initialise AR images
VAR(1:PPM.AR_P)        = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:PPM.AR_P
	VAR(i).fname   = sprintf('AR_%04d.img',i);
	VAR(i).descrip = sprintf('Autoregressive coefficient (%04d)',i);
	spm_unlink(VAR(i).fname)
end
VAR   = spm_create_vol(VAR,'noopen');

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...initialised');        %-#

% Set up masking details
%-If xM is not a structure then assumme it's a vector of thresholds
%-----------------------------------------------------------------------
try
	xM = SPM.xM;
catch
	xM = -ones(nScan,1)/0;
end

% Analyse volume/masked volume/slices ?
spm_input('Data to analyse...',1,'d')
Ctype = {
    'Volume',...
        'Masked Volume',...
        'Slices'};
str   = 'Select space';
Sel   = spm_input(str,2,'m',Ctype);
space_type = Ctype{Sel};

switch space_type,
case 'Masked Volume',
    PPM.M_X=1;
    PPM.AN_slices=[1:1:zdim];
case 'Volume',
    PPM.M_X=0;
    PPM.AN_slices=[1:1:zdim];
case 'Slices',
    PPM.M_X=0;
    PPM.AN_slices = spm_input(['Enter slice numbers eg. 3 14 2'],'+1');
end


% Explicit masking ?
%M_X = spm_input('explicitly mask images?','+1','y/n',[1,0],2);

if PPM.M_X, 
    M_P = spm_get(Inf,'*.img',{'select mask images'}); 
else, 
    M_P = {}; 
end

if ~isempty(M_P)
	xM.VM  = spm_vol(char(M_P));
end

if ~isstruct(xM)
	xM = struct(	'T',	[],...
			'TH',	xM,...
			'I',	0,...
			'VM',	{[]},...
			'xs',	struct('Masking','analysis threshold'));
end

%-Intialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM    = struct(		'fname',	'mask.img',...
			'dim',		[DIM',spm_type('uint8')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM    = spm_create_vol(VM);


%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Cycle to avoid memory problems (plane by plane)
%=======================================================================
spm_progress_bar('Init',100,'VB estimation','');
spm('Pointer','Watch')

%-maxMem is the maximum amount of data processed at a time (bytes)
%-----------------------------------------------------------------------
global defaults
MAXMEM = defaults.stats.maxmem*10;
blksz  = ceil(MAXMEM/8/nScan);
nbch     = ceil(xdim*ydim/blksz);				%-# blocks

SHp=0; % Sum of noise variance hyperparameters
xords = [1:xdim]'*ones(1,ydim); xords = xords(:)';  % plane X coordinates
yords = ones(xdim,1)*[1:ydim];  yords = yords(:)';  % plane Y coordinates
S     = 0;                                          % Number of in-mask voxels
s     = 0;                                          % Volume (voxels > UF)

% Initialise aspects of slice variables common to all slices
slice_template = spm_vb_init_volume (SPM.xX.X,PPM.AR_P);

index=1;
for  z = 1:zdim,

    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords   = z*ones(xdim*ydim,1)';	%-plane Z coordinates
    CrBl    = [];	%-conditional parameter estimates
    CrHp    = [];	% VB hyperparameter estimates
    CrPsd    = [];	% 
    CrAR    = [];	% 
    Q       = [];	%-in mask indices for this plane
    
    analyse_slice = length(find((PPM.AN_slices-z)==0));
    if analyse_slice
        for bch = 1:nbch			%-loop over bunches of lines (planks)
            
            %-# Print progress information in command window
            %---------------------------------------------------------------
            str   = sprintf('Plane %3d/%-3d, block %3d/%-3d',z,zdim,bch,nbch);
            fprintf('\r%-40s: %30s',str,' ')                                %-#
            
            %-construct list of voxels in this block
            %---------------------------------------------------------------
            I     = [1:blksz] + (bch - 1)*blksz;
            I     = I(I <= xdim*ydim);			%-truncate
            xyz   = [xords(I); yords(I); zords(I)];		%-voxel coordinates
            nVox  = size(xyz,2);
            
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
                j      = xM.VM(i).mat\M*[xyz;ones(1,nVox)];
                
                %-Load mask image within current mask & update mask
                %-------------------------------------------------------
                Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm)) > 0;
            end
            
            %-Get the data in mask, compute threshold & implicit masks
            %---------------------------------------------------------------
            Y     = zeros(nScan,nVox);
            for i = 1:nScan
                
                %-Load data in mask
                %-------------------------------------------------------
                if ~any(Cm), break, end			%-Break if empty mask
                Y(i,Cm)  = spm_get_data(VY(i),xyz(:,Cm));
                
                Cm(Cm)   = Y(i,Cm) > xM.TH(i);		%-Threshold (& NaN) mask
                if xM.I & ~YNaNrep & xM.TH(i) < 0	%-Use implicit mask
                    Cm(Cm) = abs(Y(i,Cm)) > eps;
                end
            end
            
            %-Mask out voxels where data is constant
            %---------------------------------------------------------------
            Cm(Cm) = any(diff(Y(:,Cm),1));
            Y      = Y(:,Cm);				%-Data within mask
            CrS    = sum(Cm);				%- Number of current voxels
            
            %-Conditional estimates (per partition, per voxel)
            %---------------------------------------------------------------
            beta  = zeros(nBeta,CrS);
            Hp    = zeros(1,  CrS);
            AR    = zeros(PPM.AR_P, CrS);
            Psd   = zeros(nPsd, CrS);
            
            if CrS
                
                % - Assume single session 
                vxyz  = spm_vb_neighbors(xyz(:,Cm)');
                
                slice = slice_template;
                slice.verbose=1;
                slice.update_w=1;
                slice.update_lambda=1;
                slice.update_F=0;
                
                slice=spm_vb_set_priors(slice,PPM.priors,vxyz);
                slice=spm_vb_glmar(Y,slice);
                
                if PPM.AR_P > 0
                    AR(1:PPM.AR_P,:)=slice.ap_mean;
                end
                
                if PPM.compute_F
                    slice.F=spm_vb_F(Y,slice);
                    PPM.slice(z).F=slice.F;
                end
                
                beta           = slice.wk_mean;
                Hp(1,:)        = sqrt(1./slice.mean_lambda');
                Psd            = slice.w_dev;
                
                % Get slice-wise Taylor approximation to posterior correlation
                slice = spm_vb_taylor_R (Y,slice);
                PPM.slice(z).mean=slice.mean;
                PPM.slice(z).elapsed_seconds=slice.elapsed_seconds;
                
                
                
            end % if CrS
            
            
            %-Append new inmask voxel locations and volumes
            %---------------------------------------------------------------
            XYZ(:,S + [1:CrS]) = xyz(:,Cm);		%-InMask XYZ voxel coords
            Q                  = [Q I(Cm)];		%-InMask XYZ voxel indices
            S                  = S + CrS;		%-Volume analysed (voxels)
            
            % Sum of noise variance hyperparameters
            SHp       = SHp + sum(Hp(1,:));
            
            %-Save for current plane in memory as we go along
            %---------------------------------------------------------------
            CrBl = [CrBl beta];
            CrHp = [CrHp Hp];
            CrPsd = [CrPsd Psd];
            CrAR = [CrAR AR];
            
        end % (bch)
    end % (if analyse_slice)

    %-Write Mask image
	%-------------------------------------------------------------------
	j   = sparse(xdim,ydim);
	if length(Q), j(Q) = 1; end
	VM    = spm_write_plane(VM, j, z);

    %-Write conditional beta images
    %-------------------------------------------------------------------
    j   = NaN*ones(xdim,ydim);
    for i = 1:nBeta
        if length(Q), j(Q) = CrBl(i,:); end
	    Vbeta(i)  = spm_write_plane(Vbeta(i),j,z);
    end

    %-Write SD error images
    %-------------------------------------------------------------------
    j   = NaN*ones(xdim,ydim);
	if length(Q), j(Q) = CrHp(1,:); end
	VHp    = spm_write_plane(VHp,j,z);
    
    
    %-Write posterior standard-deviation of beta images
    %-------------------------------------------------------------------
    j   = NaN*ones(xdim,ydim);
    for i = 1:nPsd
	    if length(Q), j(Q) = CrPsd(i,:); end
	    VPsd(i)    = spm_write_plane(VPsd(i),j,z);
    end
    
    %-Write AR images
    %-------------------------------------------------------------------
    j   = NaN*ones(xdim,ydim);
    for i = 1:PPM.AR_P
	    if length(Q), j(Q) = CrAR(i,:); end
	    VAR(i)    = spm_write_plane(VAR(i),j,z);
    end
    
    %-Report progress
    %-------------------------------------------------------------------
    spm_progress_bar('Set',100*(z - 1)/zdim);


end % (for z = 1:zdim)

fprintf('\n')                                                        %-#
spm_progress_bar('Clear')
spm('Pointer','Arrow');
%=======================================================================
% - P O S T   E S T I M A T I O N
%=======================================================================

if S == 0, warning('No inmask voxels - empty analysis!'), end


 %-"close" written image files, updating scalefactor information
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...closing files')      %-#
Vbeta      = spm_close_vol(Vbeta);
VHp        = spm_close_vol(VHp);
VPsd        = spm_close_vol(VPsd);
VAR        = spm_close_vol(VAR);
VM         = spm_close_vol(VM);

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#

%-Create 1st contrast for 'effects of interest' (all if not specified)
%=======================================================================
Fcname          = 'effects of interest';
try
	iX0     = [xX.iB xX.iG];
catch
	iX0     = [];
end

xX.xKXs = spm_sp('Set',spm_filter(xX.K,xX.X));		% ** Not Whitened **
xX.erdf = size(xX.X,1); % Just set to number of scans so, when 
                        % we assess the results, spm_getSPM is happy
xX.W= eye(size(xX.X,1)); % Set whitening matrix to identity -
                         % we must set it to keep spm_graph happy
                         
xCon            = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);

%-Compute scaled design matrix for display purposes
%-----------------------------------------------------------------------
xX.nKX        = spm_DesMtx('sca',xX.xKXs.X,xX.name);

%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#



%-place fields in SPM
%-----------------------------------------------------------------------
SPM.xVol.XYZ   = XYZ(:,1:S);	%-InMask XYZ coords (voxels)
SPM.xVol.M     = M;				%-voxels -> mm
SPM.xVol.iM    = inv(M);		%-mm -> voxels
SPM.xVol.DIM   = DIM;			%-image dimensions
SPM.xVol.S     = S;				%-Volume (voxels)
SPM.xVol.R     = 100;           % Set R - number of RESELS - to arbitrary value
                                % as, if R not set, SPM will think model has not 
                                % been estimated
SPM.xVol.FWHM = 10;             % Set to arbitrary value so spm_getSPM is happy
                                
SPM.VCbeta = Vbeta;			    % Filenames - parameters
SPM.VHp    = VHp;			    % Filenames - hyperparameters
SPM.VPsd   = VPsd;			    % Filenames - hyperparameters
SPM.VAR    = VAR;			    % Filenames - hyperparameters
SPM.VM     = VM;				%-Filehandle - Mask

%PPM.Cb = 1;                     % Default threshold for effect size (1 per cent)

PPM.Gamma  = 1;                 % Default threshold for effect size (1 per cent)
SPM.PPM    = PPM;			    % PPM structure

SPM.xX     = xX;                %-design structure
SPM.xM     = xM;				%-mask structure

SPM.xCon   = xCon;				%-contrast structure

save SPM SPM

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...use the results section for assessment\n\n')             %-#


