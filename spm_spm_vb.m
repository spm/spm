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
% The space to be analysed is a 'Volume', 'Slices' or 'Region'
% For 'Slices' the numbers of the slices to be analysed are then entered
% For 'Region' the centre and radius of a circular region must be specified
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
%                   .update_F, whether model evidence was computed
%
%                   .Gamma, default effect size threshold (used in spm_getSPM)
%
%                   .M_X, set to 1 if an explicit mask is used (otherwise 0)
%
%                   The following parameters are set if the space
%                   to be analysed chosen as 'Slices'
%                       .AN_slices, numbers of slices analysed
%                       .Sess(s).slice(z), further info about GLM-AR model at slice z
%                                  eg. slice(z).F is evidence for slice z 
%                                  (if computed)
%                       where s is the session number
%
%   For each session the following fields are also specified:
%
%	SPM.PPM.Sess(s).VHp        - Handle of standard deviation of the error (SDerror)
%   SPM.PPM.Sess(s).VAR        - Handles of AR coefficient images (AR_????)
%
%   If contrasts have been specified the following fields are also updated:
%
%   SPM.xCon(ic).Vcon    
%   SPM.PPM.Vcon_sd(ic)
%
%   where ic is the contrast index. 
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
% Sess%s%_SDerror.{img,hdr} 
% This is a 16-bit (float) image of the standard deviation of the error
% for session s.
% Voxels outside the analysis mask (mask.img) are given value NaN.
%
% Sess%s%_AR_????.{img,hdr} 
% This is a 16-bit (float) image of AR coefficients for session s.
% The image files are numbered according to the order of the corresponding
% AR coefficient.
% Voxels outside the analysis mask (mask.img) are given value NaN.
%
%                           ----------------
%
% Will Penny $Id$

% Let later functions know (eg. spm_contrasts) that 
% estimation was with VB
SPM.PPM.VB=1;

% Get number of sessions
nsess=length(SPM.Sess);

try 
    SPM.PPM.window;
catch
    SPM.PPM.window=1;
end

if SPM.PPM.window
    %-Say hello
    %-----------------------------------------------------------------------
    Finter   = spm('FigName','Stats: Bayesian Specification ...'); 
    spm('Pointer','Arrow');
end

if nargin ==0
    %-Get SPM.mat 
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

if SPM.PPM.window
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
			'dim',		DIM',...
			'dt',		[spm_type('float32'), spm_platform('bigend')],...
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
try 
    SPM.PPM.AR_P;
catch
    str=' AR model order';
    SPM.PPM.AR_P      = spm_input(str,1,'e',0,[1 1]);
end

% Specify type of prior
try
    SPM.PPM.priors.WA;
catch
    spm_input('Overall Priors ...',1,'d');
    Ctype = {
        'Spatial - GMRF',...
            'Spatial - LORETA',...
            'Voxel - Shrinkage',...
            'Voxel - Uninformative'};
    str   = 'Select prior';
    Sel   = spm_input(str,2,'m',Ctype);
    SPM.PPM.priors.WA = Ctype{Sel};
end

% Override option for prior on A - none
try
    SPM.PPM.priors.overA;
catch
    SPM.PPM.priors.overA='none';
end

% Compute evidence at each iteration ?
try
    SPM.PPM.update_F;
catch
    SPM.PPM.update_F=0;
end

% Initialise Error SD image(s)
for s=1:nsess,
    SPM.PPM.Sess(s).VHp        = deal(struct(...
        'fname',	[],...
        'dim',		DIM',...
	'dt,		[spm_type('float64'), spm_platform('bigend')],...
        'mat',		M,...
        'pinfo',	[1 0 0]',...
        'descrip',	''));
    
    SPM.PPM.Sess(s).VHp.fname   = sprintf('Sess%d_SDerror.img',s);
    SPM.PPM.Sess(s).VHp.descrip = sprintf('Sess%d Error SD',s);
    spm_unlink(SPM.PPM.Sess(s).VHp.fname);
    SPM.PPM.Sess(s).VHp   = spm_create_vol(SPM.PPM.Sess(s).VHp,'noopen');
end

% Initialise Posterior SD images
nPsd=size(SPM.xX.X,2);
VPsd(1:nPsd)        = deal(struct(...
			'fname',	[],...
			'dim',		DIM',...
			'dt',		[spm_type('float64'), spm_platform('bigend')],...
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
for s=1:nsess,
    SPM.PPM.Sess(s).VAR(1:SPM.PPM.AR_P)        = deal(struct(...
        'fname',	[],...
        'dim',		DIM',...
	'dt',		[spm_type('float64'), spm_platform('bigend')],...
        'mat',		M,...
        'pinfo',	[1 0 0]',...
        'descrip',	''));
    for i = 1:SPM.PPM.AR_P
        SPM.PPM.Sess(s).VAR(i).fname   = sprintf('Sess%d_AR_%04d.img',s,i);
        SPM.PPM.Sess(s).VAR(i).descrip = sprintf('Sess%d Autoregressive coefficient (%04d)',s,i);
        spm_unlink(SPM.PPM.Sess(s).VAR(i).fname)
    end
    SPM.PPM.Sess(s).VAR   = spm_create_vol(SPM.PPM.Sess(s).VAR,'noopen');
end

% Find number of pre-specified contrasts
if isfield(SPM,'xCon')
    ncon=length(SPM.xCon);
else
    ncon=0;
end

% Initialise contrast and contrast SD images
%-----------------------------------------------------------
for ic=1:ncon,
    SPM.xCon(ic).Vcon = struct(...
        'fname',  sprintf('con_%04d.img',ic),...
        'dim',    DIM',...
	'dt',     [16, spm_platform('bigend')],...
        'mat',    M,...
        'pinfo',  [1,0,0]',...
        'descrip',sprintf('SPM contrast - %d: %s',ic,SPM.xCon(ic).name));
    
    V= struct(...
        'fname',  sprintf('con_sd_%04d.img',ic),...
        'dim',    DIM',...
	'dt',     [16, spm_platform('bigend')],...
        'mat',    M,...
        'pinfo',  [1,0,0]',...
        'descrip',sprintf('PPM contrast SD - %d: %s',ic,SPM.xCon(ic).name));
    
    SPM.xCon(ic).Vcon = spm_create_vol(SPM.xCon(ic).Vcon,'noopen');
    V=spm_create_vol(V,'noopen');
    SPM.PPM.Vcon_sd(ic) = V;
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');        %-#

% Set up masking details
%-If xM is not a structure then assumme its a vector of thresholds
%-----------------------------------------------------------------------
try
	xM = SPM.xM;
catch
	xM = -ones(nScan,1)/0;
end

% Explicit masking ?
try
    SPM.PPM.M_X;
catch
    SPM.PPM.M_X = spm_input('Explicitly mask images?','+1','y/n',[1,0],2);
end

if SPM.PPM.M_X, 
    M_P = spm_get(Inf,'*.img',{'Select mask images'}); 
else, 
    M_P = {}; 
end

% Analyse volume/slices ?
try 
    SPM.PPM.space_type;
catch
    spm_input('Data to analyse...',1,'d')
    Ctype = {
        'Volume',...
            'Slices',...
            'Region'};
    str   = 'Select space';
    Sel   = spm_input(str,2,'m',Ctype);
    SPM.PPM.space_type = Ctype{Sel};
end

switch SPM.PPM.space_type,
case 'Volume',
    SPM.PPM.AN_slices=[1:1:zdim];
case 'Slices',
    try
        SPM.PPM.AN_slices;
    catch
        SPM.PPM.AN_slices = spm_input(['Enter slice numbers eg. 3 14 2'],'+1');
    end
case 'Region',
    try 
        SPM.PPM.centre;
    catch
        SPM.PPM.centre = spm_input(['Enter centre co-ordinates eg. 3 14 2'],'+1');
    end
    SPM.PPM.AN_slices=SPM.PPM.centre(3);
    try 
        SPM.PPM.radius;
    catch
        SPM.PPM.radius = spm_input(['Enter radius (in voxels) eg. 5'],'+1');
    end
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
			'dim',		DIM',...
			'dt',		[spm_type('uint8'), spm_platform('bigend')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM    = spm_create_vol(VM);


%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Cycle to avoid memory problems (plane by plane)

if SPM.PPM.window
    %=======================================================================
    spm_progress_bar('Init',100,'VB estimation','');
    spm('Pointer','Watch')
end

global defaults

SHp=0; % Sum of noise variance hyperparameters
xords = [1:xdim]'*ones(1,ydim); xords = xords(:)';  % plane X coordinates
yords = ones(xdim,1)*[1:ydim];  yords = yords(:)';  % plane Y coordinates
S     = 0;                                          % Number of in-mask voxels
s     = 0;                                          % Volume (voxels > UF)


% Initialise aspects of slice variables common to all slices
if nsess > 1
    for s=1:nsess,
        X=SPM.xX.X(SPM.Sess(s).row,SPM.Sess(s).col);
        X=[X ones(length(SPM.Sess(s).row),1)]; % Add on constant 
        slice_template(s) = spm_vb_init_volume (X,SPM.PPM.AR_P);
    end
else
    slice_template(1) = spm_vb_init_volume (SPM.xX.X,SPM.PPM.AR_P);
end

% Get matrices that will remove low-frequency drifts 
% if high pass filters have been specified
for s=1:nsess,
    sess_nScan=length(SPM.xX.K(s).row);
    if size(SPM.xX.K(s).X0,2) > 0
        X0=SPM.xX.K(s).X0;
        hpf(s).R0=eye(sess_nScan)-X0*pinv(X0);
    else
        hpf(s).R0=eye(sess_nScan);
    end    
end

% Set maximum number of VB iterations per slice
try
    SPM.PPM.maxits;
catch
    SPM.PPM.maxits=4;
end
for s=1:nsess,
    slice_template(s).maxits=SPM.PPM.maxits;
    slice_template(s).verbose=1;
    slice_template(s).update_w=1;
    slice_template(s).update_lambda=1;
    slice_template(s).update_F=SPM.PPM.update_F;
end

index=1;
for  z = 1:zdim,

    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords   = z*ones(xdim*ydim,1)';	%-plane Z coordinates
    CrBl    = [];	%-conditional parameter estimates
    CrPsd    = [];	% 
    for s=1:nsess,
        Sess(s).CrAR    = [];	% AR estimates
        Sess(s).CrHp    = [];	% 
    end
    Cr_con  = [];   % Contrasts
    Cr_con_var = []; % Contrast variances
    Q       = [];	%-in mask indices for this plane
    
    analyse_slice = length(find((SPM.PPM.AN_slices-z)==0));
    if analyse_slice
        % Only ever use 1 bunch of lines (planks)
        bch=1; nbch=1;
        
        %-# Print progress information in command window
        %---------------------------------------------------------------
        str   = sprintf('Plane %3d/%-3d',z,zdim);
        fprintf('\r%-40s: %30s',str,' ')                                %-#
        
        %-construct list of voxels 
        %---------------------------------------------------------------
        I     = [1:xdim*ydim];
        xyz   = [xords(I); yords(I); zords(I)];		%-voxel coordinates
        nVox  = size(xyz,2);
        
        %-Get data & construct analysis mask
        %===============================================================
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read & mask data')%-#
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
        CrS    = sum(Cm);				%- Number of current voxels
        
        if CrS
            
            if strcmp(SPM.PPM.space_type,'Region')
                % Update current mask to specify only this region
                d=[xyz(1,:)-SPM.PPM.centre(1);
                    xyz(2,:)-SPM.PPM.centre(2)];
                Cm(find(sum(d.^2) > SPM.PPM.radius))=0;
            end
            
            CrS = sum(Cm);
            Y   = Y(:,Cm);				%-Data within mask
            
            %-Conditional estimates (per partition, per voxel)
            %---------------------------------------------------------------
            beta  = zeros(nBeta,CrS);
            Psd   = zeros(nPsd, CrS);
            
            for s=1:nsess,
                Sess(s).Hp    = zeros(1,  CrS);
                Sess(s).AR    = zeros(SPM.PPM.AR_P, CrS);
            end
            
            if ncon > 0
                con   = zeros(ncon,CrS);
                con_var = zeros(ncon,CrS);
            end
            vxyz  = spm_vb_neighbors(xyz(:,Cm)');
            
            % Estimate model for each session separately
            for s=1:nsess;
                
                slice = slice_template(s);
                slice = spm_vb_set_priors(slice,SPM.PPM.priors,vxyz);
                
                % Filter data to remove low frequencies
                R0Y=hpf(s).R0*Y(SPM.Sess(s).row,:);
                
                % Fit model
                slice = spm_vb_glmar(R0Y,slice);
                
                % Report AR values 
                if SPM.PPM.AR_P > 0
                    % Report AR values averaged over sessions
                    %AR(1:SPM.PPM.AR_P,:)=AR(1:SPM.PPM.AR_P,:) + slice.ap_mean/nsess;
                    
                    % Report session specific AR values
                    Sess(s).AR(1:SPM.PPM.AR_P,:)=slice.ap_mean;
                end
                
                if SPM.PPM.update_F
                    SPM.PPM.Sess(s).slice(z).F=slice.F;
                end
                
                % Update regression coefficients 
                ncols=length(SPM.Sess(s).col);
                beta(SPM.Sess(s).col,:) = slice.wk_mean(1:ncols,:);
                if ncols==0
                    % Design matrix empty except for constant
                    mean_col_index=s;
                else
                    mean_col_index=SPM.Sess(nsess).col(end)+s;
                end
                beta(mean_col_index,:) = slice.wk_mean(ncols+1,:); % Session mean
                
                % Report noise variances averaged over sessions
                % Hp(1,:)        = Hp(1,:)+(1/nsess)*sqrt(1./slice.mean_lambda');
                
                % Report session-specific noise variances 
                Sess(s).Hp(1,:)        = sqrt(1./slice.mean_lambda');
                
                % Store regression coefficient posterior standard deviations 
                Psd (SPM.Sess(s).col,:) = slice.w_dev(1:ncols,:);
                Psd (mean_col_index,:) = slice.w_dev(ncols+1,:);
                
                % Update contrast variances
                if ncon > 0
                    for ic=1:ncon,
                        CC=SPM.xCon(ic).c;
                        % Get relevant columns of contrast
                        CC=[CC(SPM.Sess(s).col) ; 0];
                        for i=1:CrS,
                            con_var(ic,i)=con_var(ic,i)+CC'*slice.w_cov{i}*CC;
                        end
                    end
                end
                
                % Get slice-wise Taylor approximation to posterior correlation
                slice = spm_vb_taylor_R (R0Y,slice);
                SPM.PPM.Sess(s).slice(z).mean=slice.mean;
                SPM.PPM.Sess(s).slice(z).elapsed_seconds=slice.elapsed_seconds;
                
                % Save Coefficient RESELS and number of voxels
                SPM.PPM.Sess(s).slice(z).gamma_tot=slice.gamma_tot;
                SPM.PPM.Sess(s).slice(z).N=slice.N;
            
                clear slice;
            end % loop over sessions
            
            % Get contrasts 
            if ncon > 0
                for ic=1:ncon,
                    CC=SPM.xCon(ic).c;
                    con(ic,:)=CC'*beta;
                end
            end
            
            %-Append new inmask voxel locations and volumes
            %---------------------------------------------------------------
            XYZ(:,S + [1:CrS]) = xyz(:,Cm);		%-InMask XYZ voxel coords
            Q                  = [Q I(Cm)];		%-InMask XYZ voxel indices
            S                  = S + CrS;		%-Volume analysed (voxels)
            
            % Sum of noise variance hyperparameters
            SHp       = SHp + (1/nsess)*sum(Sess(s).Hp(1,:));
            
            %-Save for current plane in memory as we go along
            %---------------------------------------------------------------
            CrBl = [CrBl beta];
            CrPsd = [CrPsd Psd];
            for s=1:nsess,
                Sess(s).CrHp = [Sess(s).CrHp Sess(s).Hp];
                Sess(s).CrAR = [Sess(s).CrAR Sess(s).AR];
            end
            if ncon > 0
                Cr_con = [Cr_con con];
                Cr_con_var = [Cr_con_var con_var];
            end
        
        end % if CrS
        
    end % (if analyse_slice)

    %-Write Mask image
	%-------------------------------------------------------------------
	j   = sparse(xdim,ydim);
	if length(Q), j(Q) = 1; end
	VM    = spm_write_plane(VM, j, z);

    %-Write conditional beta images
    %-------------------------------------------------------------------
    j   = repmat(NaN,xdim,ydim);
    for i = 1:nBeta
        if length(Q), j(Q) = CrBl(i,:); end
	    Vbeta(i)  = spm_write_plane(Vbeta(i),j,z);
    end

    %-Write SD error images
    %-------------------------------------------------------------------
    for s=1:nsess,
        j = repmat(NaN,xdim,ydim);
        if length(Q), j(Q) = Sess(s).CrHp(1,:); end
        SPM.PPM.Sess(s).VHp    = spm_write_plane(SPM.PPM.Sess(s).VHp,j,z);
    end
    
    %-Write posterior standard-deviation of beta images
    %-------------------------------------------------------------------
    j   = repmat(NaN,xdim,ydim);
    for i = 1:nPsd
	    if length(Q), j(Q) = CrPsd(i,:); end
	    VPsd(i)    = spm_write_plane(VPsd(i),j,z);
    end
    
    %-Write AR images
    %-------------------------------------------------------------------
    for s=1:nsess,
        j   = repmat(NaN,xdim,ydim);
        for i = 1:SPM.PPM.AR_P
            if length(Q), j(Q) = Sess(s).CrAR(i,:); end
            SPM.PPM.Sess(s).VAR(i)    = spm_write_plane(SPM.PPM.Sess(s).VAR(i),j,z);
        end
    end
    
    
    % Write contrast and contrast SD images
    if ncon > 0
        j   = repmat(NaN,xdim,ydim);
        for ic=1:ncon
            if length(Q), j(Q) = Cr_con(ic,:); end
	        SPM.xCon(ic).Vcon    = spm_write_plane(SPM.xCon(ic).Vcon,j,z);
        end
        j   = repmat(NaN,xdim,ydim);
        for ic=1:ncon
            if length(Q), j(Q) = sqrt(Cr_con_var(ic,:)); end
	        SPM.PPM.Vcon_sd(ic)    = spm_write_plane(SPM.PPM.Vcon_sd(ic),j,z);
        end
    end

    if SPM.PPM.window
        %-Report progress
        %-------------------------------------------------------------------
        spm_progress_bar('Set',100*(z - 1)/zdim);
    end

end % (for z = 1:zdim)

fprintf('\n')   %-#
if SPM.PPM.window
    spm_progress_bar('Clear')
    spm('Pointer','Arrow');
end

%=======================================================================
% - P O S T   E S T I M A T I O N
%=======================================================================

if S == 0, warning('No inmask voxels - empty analysis!'), end


 %-"close" written image files, updating scalefactor information
%=======================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...closing files')      %-#
Vbeta      = spm_close_vol(Vbeta);
VM         = spm_close_vol(VM);
VPsd        = spm_close_vol(VPsd);
for s=1:nsess,
    SPM.PPM.Sess(s).VAR        = spm_close_vol(SPM.PPM.Sess(s).VAR);
    SPM.PPM.Sess(s).VHp        = spm_close_vol(SPM.PPM.Sess(s).VHp);
end
% Close contrast and contrast variance images
if ncon > 0
    for ic=1:ncon,
        SPM.xCon(ic).Vcon    = spm_close_vol(SPM.xCon(ic).Vcon);
        SPM.PPM.Vcon_sd(ic)    = spm_close_vol(SPM.PPM.Vcon_sd(ic));
    end
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')               %-#

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
SPM.VPsd   = VPsd;			    % Filenames - hyperparameters
SPM.VM     = VM;				%-Filehandle - Mask

SPM.PPM.Gamma  = 1;             % Default threshold for effect size (1 per cent)

SPM.xX     = xX;                %-design structure
SPM.xM     = xM;				%-mask structure

save SPM SPM

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')               %-#


if SPM.PPM.window
    %=======================================================================
    %- E N D: Cleanup GUI
    %=======================================================================
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
    fprintf('...use the results section for assessment\n\n')             %-#
end


