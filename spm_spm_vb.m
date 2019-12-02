  function [SPM] = spm_spm_vb(SPM)
% VB estimation of GLM-AR models with spatial regularisation
% FORMAT [SPM] = spm_spm_vb(SPM)
%
% This function implements a VB estimation scheme for GLM-AR models. 
% Both regression coefficients and AR coefficients are spatially
% regularised. The algorithm is described in a series of papers:
%
% Paper VB1: W. Penny, S. Kiebel and K. Friston (2003).
%            Variational Bayesian inference for fMRI time series. 
%            NeuroImage 19(3), pp 727-741.
%
% Paper VB2: W. Penny, N. Trujillo-Barreto and K. Friston (2005).
%            Bayesian fMRI time series analysis with spatial priors.
%            NeuroImage 24(2), pp 350-362.
%
% Paper VB3: W. Penny and G. Flandin (2005).
%            Bayesian analysis of single-subject fMRI: SPM implementation.
%            Technical Report. WDIN, UCL.
%
% Paper VB4: W. Penny, G. Flandin and N. Trujillo-Barreto (2007).
%            Bayesian comparison of spatially regularised general linear
%            models. Human Brain Mapping 28(4):275-293.
%
% Paper VB5: W. Penny and G. Flandin (2005).
%            Bayesian analysis of fMRI data with spatial priors.
%            2005 Proceedings of the Joint Statistical Meeting, Section on
%            Statistical Graphics [CDROM], Alexandria, VA: American 
%            Statistical Association.  
%
% The space to be analysed is a 'Volume', 'Slices' or 'Clusters'.
% For 'Slices' the numbers of the slices to be analysed are then entered.
% For 'Clusters' a mask specifying the volume to be analysed is entered.
% ______________________________________________________________________ 
%
% Required fields of SPM:
%
% xY.VY - nScan x 1 struct array of mapped image volumes
%         Images must have the same orientation, voxel size and data type
%       - Any scaling should have already been applied via the image
%         handle scalefactors.
%
% xX    - Structure containing design matrix information
%       - Required fields are:
%         xX.X      - Design matrix (raw, not temporally smoothed)
%         xX.name   - cellstr of parameter names corresponding to 
%                     columns of design matrix
%
% xM    - Structure containing masking information, or a simple column
%         vector of thresholds corresponding to the images in VY.
%       - If a structure, the required fields are:
%         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1 => implicit zero/NaN mask)
%         xM.VM - struct array of mapped explicit mask image volumes
%       - (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or 
%                 data type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%       - Note that voxels with constant data (i.e. the same value 
%         across scans) are also automatically masked out.
%         
% ______________________________________________________________________
%
% spm_spm_vb adds the following fields to SPM:
%
%   SPM.VCbeta - Handles of posterior parameter estimates 
%                (Cbeta_????)
%   SPM.VPsd   - Handles of SD of posterior parameter estimates
%                (SDbeta_????)
%
%   SPM.PPM    - Posterior Probability Map data structure
%
%                .VB=1, tells later functions (spm_contrasts, spm_graph)
%                       that parameters were estimated using VB
%                       
%                .AR_P, assumed AR model order
%
%                .priors, type of priors used (eg. 'Spatial-GMRF')
%                       see spm_vb_set_priors.m
%
%                .update_F, whether model evidence is to be computed
%
%                .Gamma, default effect size threshold (used in spm_getSPM)
%
%                   info is stored for each "block", where a block
%                   is either a slice or subvolume, computed using a 
%                   graph-partitioning algorithm. This is stored in
%                       .Sess(s).block(z), further info about GLM-AR 
%                                  model at block z eg. block(z).F is
%                                  evidence for block z (if computed)
%                                  where s is the session number
%                   The following parameters are set if the space to be
%                   analysed chosen as 'Slices'
%                       .AN_slices, numbers of slices analysed
%
% For each session the following fields are also specified:
%
%   SPM.PPM.Sess(s).VHp  - Handle of standard deviation of the error
%                          (Sess%s%_SDerror)
%   SPM.PPM.Sess(s).VAR  - Handles of AR coefficient images
%                          (Sess%s%_AR_????)
%
% If contrasts have been specified prior to estimation (this is 
% recommended) the following fields are also updated:
%
%   SPM.xCon(ic).Vcon    
%   SPM.PPM.Vcon_sd(ic)
%
% where ic is the contrast index. 
% ______________________________________________________________________
%    
%
% The following images are written to file:
%
% mask.{img,hdr}                                - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
% Note mask.<ext> is only written if the selected space is 'Volume' or 
% 'Masked Volume' (ie not 'Slices')
%
% labels.<ext>                                - block labels
% 8-bit (uint8) image of zero-s & integers from 1 to max no. of blocks, 
% e.g. slices or subvoumes, indicating which block a voxel belongs. 
% This info is also stored in SPM.xVol.labels (same order as XYZ matrix),  
% for all analysis space options.
%
% Cbeta_????.<ext>  
% These are 16-bit (float) images of the parameter posteriors. The image
% files are numbered according to the corresponding column of the
% design matrix.
% Voxels outside the analysis mask (mask.<ext>) are given value NaN.
%
% SDbeta_????.<ext> 
% These are 16-bit (float) images of the standard deviation of parameter
% posteriors. 
% The image files are numbered according to the corresponding column of
% the design matrix.
% Voxels outside the analysis mask (mask.<ext>) are given value NaN.
%
% Sess%s%_SDerror.<ext>
% This is a 16-bit (float) image of the standard deviation of the error
% for session s.
% Voxels outside the analysis mask (mask.<ext>) are given value NaN.
%
% Sess%s%_AR_????.<ext>
% These are 16-bit (float) images of AR coefficients for session s.
% The image files are numbered according to the order of the 
% corresponding AR coefficient.
% Voxels outside the analysis mask (mask.<ext>) are given value NaN.
%_______________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% Will Penny, Nelson Trujillo-Barreto and Lee Harrison
% $Id: spm_spm_vb.m 7739 2019-12-02 14:00:18Z guillaume $


%-Get SPM.mat if necessary
%-----------------------------------------------------------------------
if ~nargin
    [Pf, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd = spm_file(Pf,'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd = swd;
end

%-Change to SPM.swd if specified
%-----------------------------------------------------------------------
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

%-Let later functions (spm_contrasts, spm_graph) know that estimation
% was with Variational Bayes
%-----------------------------------------------------------------------
SPM.PPM.VB = 1;

%-Display
%-----------------------------------------------------------------------
try 
    SPM.PPM.window;
catch
    SPM.PPM.window = 1;
end

if SPM.PPM.window
    %-Say hello
    %-------------------------------------------------------------------
    Finter = spm('FigName','Stats: Bayesian Estimation ...'); 
    spm('Pointer','Arrow');
end

%-Delete files from previous analyses
%-----------------------------------------------------------------------
if SPM.PPM.window
    if ~isempty(spm_select('List',SPM.swd,'^mask\..{3}$'))
        str   = {'Current directory contains SPM estimation files:',...
                 'pwd = ',SPM.swd,...
                 'Existing results will be overwritten!'};
        
        abort = spm_input(str,1,'bd','stop|continue',[1,0],1);
        if abort
            spm('FigName','Stats: done',Finter); spm('Pointer','Arrow');
            return
        else
            str = sprintf('Overwriting old results\n\t (pwd = %s) ',SPM.swd);
            warning(str)
            drawnow
        end
    end
end

fspm = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
        '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
        '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};

fppm = {'^Cbeta_.{4}\..{3}$', '^LogEv\..{3}$', '^Sess.+_SDerror\..{3}$',...
        '^SDbeta_.{4}\..{3}$', '^Sess.+_AR_.{4}\..{3}$', '^con_sd_.{4}\..{3}$'};

files = {fspm{:} fppm{:}};

for i=1:length(files)
    j = spm_select('List',pwd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end

%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Get number of sessions
%-----------------------------------------------------------------------
nsess  = length(SPM.Sess);

%-Get image dimensions and data
%-----------------------------------------------------------------------
VY     = SPM.xY.VY;
M      = VY(1).mat;
DIM    = VY(1).dim(1:3)';
xdim   = DIM(1); ydim = DIM(2); zdim = DIM(3);

%-Get design matrix
%-----------------------------------------------------------------------
xX     = SPM.xX;
[nScan,nBeta] = size(xX.X);
nPsd   = nBeta;

%-Find number of pre-specified contrasts
%-----------------------------------------------------------------------
try
    ncon = length(SPM.xCon);
catch
    ncon = 0;
end

%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising');               %-#

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ   = zeros(3,xdim*ydim*zdim);
labels = zeros(1,xdim*ydim*zdim);

%-Initialise conditional estimate image files
%-----------------------------------------------------------------------
Vbeta(1:nBeta) = deal(struct(...
            'fname',    '',...
            'dim',      DIM',...
            'dt',       [spm_type('float32') spm_platform('bigend')],...   
            'mat',      M,...
            'pinfo',    [1 0 0]',...
            'descrip',  ''));
        
for i = 1:nBeta
    Vbeta(i).fname   = [sprintf('Cbeta_%04d',i) spm_file_ext];
    Vbeta(i).descrip = sprintf('Posterior mean of beta (%04d) - %s',i,xX.name{i});
end
Vbeta = spm_create_vol(Vbeta);

%-Initialise Posterior SD image files
%-----------------------------------------------------------------------
VPsd(1:nPsd) = deal(struct(...
    'fname',    '',...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  ''));
    
for i = 1:nPsd
    VPsd(i).fname   = [sprintf('SDbeta_%04d',i) spm_file_ext];
    VPsd(i).descrip = sprintf('Posterior SD of beta (%04d)',i);
end
VPsd = spm_create_vol(VPsd);

%-Initialise Error SD image(s)
%-----------------------------------------------------------------------
for s = 1:nsess
    SPM.PPM.Sess(s).VHp = struct(...
        'fname',    [],...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],... 
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'descrip',  '');
    
    SPM.PPM.Sess(s).VHp.fname   = [sprintf('Sess%d_SDerror',s) spm_file_ext];
    SPM.PPM.Sess(s).VHp.descrip = sprintf('Sess%d Error SD',s);
    SPM.PPM.Sess(s).VHp = spm_create_vol(SPM.PPM.Sess(s).VHp);
end

%-Initialise hyperparameter (AR 1..p and noise variance) image files
%-----------------------------------------------------------------------

%-Set number of AR coefficients
try 
    SPM.PPM.AR_P;
catch
    SPM.PPM.AR_P = 3;
end

for s=1:nsess
    for i=1:SPM.PPM.AR_P
        SPM.PPM.Sess(s).VAR(i) = struct(...
            'fname',    '',...
            'dim',      DIM',...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      M,...
            'pinfo',    [1 0 0]',...
            'descrip',  '',...
            'n',        1,...
            'private',  []);
        SPM.PPM.Sess(s).VAR(i).fname   = ...
            [sprintf('Sess%d_AR_%04d',s,i) spm_file_ext];
        SPM.PPM.Sess(s).VAR(i).descrip = ...
            sprintf('Sess%d Autoregressive coefficient (%04d)',s,i);
    end
    if SPM.PPM.AR_P > 0
        SPM.PPM.Sess(s).VAR = spm_create_vol(SPM.PPM.Sess(s).VAR);
    end
end

%-Initialise contribution map
% (voxel-wise contribution to log evidence)
%-----------------------------------------------------------------------

%-Compute evidence at each iteration ?
try
    SPM.PPM.update_F;
catch
    SPM.PPM.update_F = 0;
end

if SPM.PPM.update_F
    SPM.PPM.LogEv = struct(...
        'fname',    ['LogEv' spm_file_ext],...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],... 
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'descrip',  'Map of contribution to log-evidence');
        
    SPM.PPM.LogEv = spm_create_vol(SPM.PPM.LogEv);
end

%-Initialise contrast and contrast SD images
%-----------------------------------------------------------------------
for ic = 1:ncon
    SPM.xCon(ic).Vcon = struct(...
        'fname',   [sprintf('con_%04d',ic) spm_file_ext],...
        'dim',     DIM',...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     M,...
        'pinfo',   [1,0,0]',...
        'descrip', sprintf('SPM contrast - %d: %s',ic,SPM.xCon(ic).name));
    
    V = struct(...
        'fname',   [sprintf('con_sd_%04d',ic) spm_file_ext],...
        'dim',     DIM',...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     M,...
        'pinfo',   [1,0,0]',...
        'descrip', sprintf('PPM contrast SD - %d: %s',ic,SPM.xCon(ic).name));
    
    SPM.xCon(ic).Vcon = spm_create_vol(SPM.xCon(ic).Vcon);
    V = spm_create_vol(V);
    SPM.PPM.Vcon_sd(ic) = V;
end

% Set up masking details
%-If xM is not a structure then assumme it's a vector of thresholds
%-----------------------------------------------------------------------
try
    xM = SPM.xM;
catch
    xM = repmat(-Inf,nScan,1);
end

if ~isstruct(xM)
    xM = struct(...
        'T',   [],...
        'TH',  xM,...
        'I',   0,...
        'VM',  {[]},...
        'xs',  struct('Masking','analysis threshold'));
end

%-Initialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM = struct(...
    'fname',    ['mask' spm_file_ext],...
    'dim',      DIM',...
    'dt',       [spm_type('uint8') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm:resultant analysis mask');
VM = spm_create_vol(VM);

%=======================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');        %-#

%-Specify type of prior for regression coefficients
%-----------------------------------------------------------------------
try
    SPM.PPM.priors.W;
catch
    SPM.PPM.priors.W = 'Spatial - UGL';
end

%-Specify type of prior for AR coefficients
%-----------------------------------------------------------------------
try
    SPM.PPM.priors.A;
catch
    SPM.PPM.priors.A = 'Spatial - UGL';
end

if  ~isdeployed && strcmp(SPM.PPM.priors.A,'Robust')
    addpath(fullfile(spm('Dir'),'toolbox','mixture'));
end

%-Get structural info if necessary
%-----------------------------------------------------------------------
if strcmp(SPM.PPM.priors.A,'Discrete')
    try
        SPM.PPM.priors.SY;
    catch
        SPM.PPM.priors.SY = spm_select([1 Inf],'image',...
            'Select structural images eg. brain or grey/white/CSF'); 
    end
    SPM.PPM.priors.Sin = size(SPM.PPM.priors.SY,1);
    for j=1:SPM.PPM.priors.Sin
        xDiscrete(j)  = spm_vol(deblank(SPM.PPM.priors.SY(j,:)));
    end
end

%-Analysis space (volume/slices/clusters)
%-----------------------------------------------------------------------
try 
    SPM.PPM.space_type;
catch
    SPM.PPM.space_type = 'Volume';
end

switch lower(SPM.PPM.space_type)
    case 'volume',
        SPM.PPM.AN_slices = [1:1:zdim];
    case 'slices',
        try
            SPM.PPM.AN_slices;
        catch
            SPM.PPM.AN_slices = spm_input(['Enter slice numbers eg. 3 14 2'],'+1');
        end
    case {'clusters'}
        %-Cluster mask 
        %-----------------------------------------------------------------------
        CM = spm_vol(SPM.PPM.clustermask{1});
        SPM.PPM.AN_slices = [1:zdim];
    otherwise
        error('Unknown analysis space.');
end

%-Initialise image containing labels of each block (slice or subvolume) 
%-----------------------------------------------------------------------
VLabel = struct(...
    'fname',    ['labels' spm_file_ext],...
    'dim',      DIM',...
    'dt',       [spm_type('uint8') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'labels used to partition a volume');
VLabel = spm_create_vol(VLabel);

[xords,yords] = ndgrid(1:xdim,1:ydim);
xords = xords(:)';  % plane X coordinates
yords = yords(:)';  % plane Y coordinates
S     = 0;          % Number of in-mask voxels
s     = 0;          % Volume (voxels > UF)

%-Initialise aspects of block variables common to all blocks
%-----------------------------------------------------------------------
if nsess > 1
    for s=1:nsess
        X = SPM.xX.X(SPM.Sess(s).row,SPM.Sess(s).col);
        X = [X ones(length(SPM.Sess(s).row),1)]; % Add on constant 
        block_template(s) = spm_vb_init_volume(X,SPM.PPM.AR_P);
    end
else
    block_template(1) = spm_vb_init_volume(SPM.xX.X,SPM.PPM.AR_P);
end

%-Get matrices that will remove low-frequency drifts 
% if high pass filters have been specified
%-----------------------------------------------------------------------
for s=1:nsess
    sess_nScan = length(SPM.xX.K(s).row);
    if size(SPM.xX.K(s).X0,2) > 0
        X0 = SPM.xX.K(s).X0;
        hpf(s).R0 = eye(sess_nScan)-X0*pinv(X0);
    else
        hpf(s).R0 = eye(sess_nScan);
    end    
end

%-Set maximum number of VB iterations per block
%-----------------------------------------------------------------------
try
    SPM.PPM.maxits;
catch
    SPM.PPM.maxits=4;
end
try
    SPM.PPM.tol;
catch
    SPM.PPM.tol=0.0001;
end
try
    SPM.PPM.compute_det_D;
catch
    SPM.PPM.compute_det_D=0;
end

for s=1:nsess
    block_template(s).maxits        = SPM.PPM.maxits;
    block_template(s).tol           = SPM.PPM.tol;
    block_template(s).compute_det_D = SPM.PPM.compute_det_D;
    block_template(s).verbose       = 0;
    block_template(s).update_w      = 1;
    block_template(s).update_lambda = 1;
    block_template(s).update_F      = SPM.PPM.update_F;
end

%-Compute mask volume - before analysis
%-----------------------------------------------------------------------
fprintf('%-40s: %30s','Calculating mask',' ')                     %-#
for z = 1:zdim

    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords = repmat(z,1,xdim*ydim); %-plane Z coordinates
    Q          = [];  %-in mask indices for this plane

    if ismember(z,SPM.PPM.AN_slices)

        %-Print progress information in command window
        %---------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',z,zdim)) %-#

        %-Construct list of voxels
        %---------------------------------------------------------------
        I     = [1:xdim*ydim];
        xyz   = [xords(I); yords(I); zords(I)];      %-voxel coordinates
        nVox  = size(xyz,2);

        %-Get data & construct analysis mask
        %---------------------------------------------------------------
        Cm    = true(1,nVox);                        %-current mask

        %-Compute explicit mask
        % (note that these may not have same orientations)
        %---------------------------------------------------------------
        for i = 1:length(xM.VM)

            %-Coordinates in mask image
            %-----------------------------------------------------------
            j      = xM.VM(i).mat\M*[xyz;ones(1,nVox)];

            %-Load mask image within current mask & update mask
            %-----------------------------------------------------------

            Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm)) > 0;

        end

        if strcmp(SPM.PPM.space_type,'clusters')
            %-Coordinates in cluster mask image
            %-----------------------------------------------------------
            j      = CM.mat\M*[xyz;ones(1,nVox)];

            %-Load mask image within current mask & update mask
            %-----------------------------------------------------------

            Cm(Cm) = spm_get_data(CM,j(:,Cm)) > 0;
        end

        %-Get the data in mask, compute threshold & implicit masks
        %---------------------------------------------------------------
        Y     = zeros(nScan,nVox);
        for i = 1:nScan

            %-Load data in mask
            %-----------------------------------------------------------
            if ~any(Cm), break, end             %-Break if empty mask
            Y(i,Cm)  = spm_get_data(VY(i),xyz(:,Cm));

            Cm(Cm)   = Y(i,Cm) > xM.TH(i);      %-Threshold (& NaN) mask
            if xM.I && xM.TH(i) < 0             %-Use implicit mask
                Cm(Cm) = abs(Y(i,Cm)) > eps;
            end
        end

        %-Mask out voxels where data is constant
        %---------------------------------------------------------------
        Cm(Cm) = any(diff(Y(:,Cm),1));

        CrS = sum(Cm);
        
        if CrS,
            %-Remove isolated nodes (mask is then the same for slice 
            % and graph-partitioned analyses)
            %-----------------------------------------------------------
            vxyz = spm_vb_neighbors(xyz(:,Cm)',0);
            if any(sum(vxyz,2)==0)
                Cm(Cm) = (sum(vxyz,2)>0);
            end
        end       
        
        %-Append new inmask voxel locations and volumes
        %---------------------------------------------------------------
        Q                  = I(Cm);   %-InMask XYZ voxel indices

    end

    %-Write Mask image
    %-------------------------------------------------------------------
    j  = sparse(xdim,ydim);
    if length(Q), j(Q) = 1; end
    VM = spm_write_plane(VM, j, z);

end
fprintf('\n')

%-Remove small clusters - removes clusters containing < 16 voxels 
%-----------------------------------------------------------------------
mask        = spm_read_vols(spm_vol(VM));
[Cl,nCl]    = spm_bwlabel(mask,6);
ncl         = histc(Cl(:),[1:max(Cl(:))])';
if any(ncl < 16)
    incl = find(ncl < 16);
    for j = 1:length(incl),
        mask(find(Cl==incl(j))) = 0;
    end
    VM          = spm_write_plane(VM, mask, ':');
    [Cl,nCl]    = spm_bwlabel(mask,6);
    ncl         = histc(Cl(:),[1:max(Cl(:))])';
end

%-Compute labels
%-----------------------------------------------------------------------
nLb = 0;
Lb = zeros(size(Cl));
if strcmp(SPM.PPM.block_type,'subvolumes') % using graph partitioning
    vol     = 1;
    CUTOFF  = 1000; % minimal number of voxels in a block
    for num = 1:nCl
        I = find(Cl==num);
        if ncl(num) > CUTOFF
            N       = ncl(num);
            [x,y,z] = ind2sub([DIM(1),DIM(2),DIM(3)],I);
            xyz     = [x,y,z];
            vxyz    = spm_vb_neighbors(xyz,1);
            [edges,weights] = spm_vb_edgeweights(vxyz);
            W       = spm_vb_adjacency(edges,weights,N);
            lbs     = zeros(N,1);
            ind     = [1:N]';
            depth   = 1;
            lbs     = spm_vb_graphcut(lbs,ind,I,W,depth,'random',CUTOFF,DIM);
            lbs     = lbs + 1;
            nl      = max(lbs);
            for z = 1:nl,
                Lb(I(lbs==z)) = nLb + z;
            end
        else
            nl      = 1;
            Lb(I)   = nLb + nl;
        end
        nLb = nLb + nl;
    end
else
    % label slices
    vol     = 0;
    mask    = spm_read_vols(spm_vol(VM));
    for z = 1:zdim,
        if any(any(mask(:,:,z)))
            nLb         = nLb + 1;
            Lb(:,:,z)   = mask(:,:,z)*nLb;
        end
    end
end
nlb     = histc(Lb(:),[1:max(Lb(:))])';

%-Write VLabel
%-----------------------------------------------------------------------
VLabel  = spm_write_plane(VLabel,Lb,':'); 
    
%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================
if SPM.PPM.window
    spm_progress_bar('Init',100,'VB estimation','');
    spm('Pointer','Watch')
end

%-Block-wise analysis (a block is either a slice or sub-volume)
%-----------------------------------------------------------------------
for z = 1:nLb

    %-Print progress information in command window
    %-------------------------------------------------------------------
    str   = sprintf('Block %3d/%-3d',z,nLb);
    fprintf('\r%-40s: %30s',str,' ')                                %-#

    %-Construct list of voxels
    %-------------------------------------------------------------------
    Q           = find(Lb==z);
    [xx,yy,zz]  = ind2sub(size(mask),Q);
    xyz         = [xx,yy,zz]';                      %-voxel coordinates
    nVox        = size(xyz,2);

    %-Get data
    %-------------------------------------------------------------------
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...read & mask data');%-#

    Y     = zeros(nScan,nVox);
    for i = 1:nScan

        Y(i,:)  = spm_get_data(VY(i),xyz);

    end

    vxyz = spm_vb_neighbors(xyz',vol);

    %-Conditional estimates (per partition, per voxel)
    %-------------------------------------------------------------------
    beta  = zeros(nBeta,nVox);
    Psd   = zeros(nPsd, nVox);
    LogEv = zeros(1,nVox);

    for s=1:nsess
        Sess(s).Hp = zeros(1, nVox);
        Sess(s).AR = zeros(SPM.PPM.AR_P, nVox);
    end

    if ncon > 0
        con     = zeros(ncon,nVox);
        con_var = zeros(ncon,nVox);
    end

    %-Get structural info for that block
    %-------------------------------------------------------------------
    if strcmp(SPM.PPM.priors.A,'Discrete')
        sxyz = xDiscrete(1).mat\M*[xyz;ones(1,nVox)];
        gamma = [];
        for j=1:SPM.PPM.priors.Sin
            gamma(:,j) = spm_get_data(xDiscrete(j),sxyz)';
        end
        if SPM.PPM.priors.Sin==1
            unaccounted=find(gamma==0);
            if length(unaccounted) > 0
                % Create extra category
                SPM.PPM.priors.S=2;
                gamma(:,2)=zeros(nVox,1);
                gamma(unaccounted,2)=1;
                gamma(:,1)=ones(nVox,1)-gamma(:,2);
                SPM.PPM.priors.gamma=gamma;
            else
                SPM.PPM.priors.S=1;
                SPM.PPM.priors.gamma=ones(nVox,1);
            end
        else
            unaccounted=find(sum(gamma')==0);
            if length(unaccounted)>0
                % Create extra category
                SPM.PPM.priors.S=SPM.PPM.priors.Sin+1;
                gamma(:,SPM.PPM.priors.S)=zeros(1,nVox);
                gamma(unaccounted,SPM.PPM.priors.S)=1;
            else
                SPM.PPM.priors.S=SPM.PPM.priors.Sin;
            end
            % convert probabilities to discrete values
            SPM.PPM.priors.gamma=zeros(nVox,SPM.PPM.priors.S);
            [yy,ii]=max(gamma');
            for j=1:SPM.PPM.priors.S
                SPM.PPM.priors.gamma(find(ii==j),j)=1;
            end
        end
    end

    %-Estimate model for each session separately
    %-------------------------------------------------------------------
    for s = 1:nsess

        fprintf('Session %d',s);                                %-#
        block = block_template(s);
        block = spm_vb_set_priors(block,SPM.PPM.priors,vxyz);

        %-Filter data to remove low frequencies
        %---------------------------------------------------------------
        R0Y = hpf(s).R0*Y(SPM.Sess(s).row,:);

        %-Fit model
        %---------------------------------------------------------------
        switch SPM.PPM.priors.A
            case 'Robust',
                %k=SPM.PPM.priors.k;
                block = spm_vb_robust(R0Y,block);
            otherwise
                block = spm_vb_glmar(R0Y,block);
        end

        %-Report AR values
        %---------------------------------------------------------------
        if SPM.PPM.AR_P > 0
            % session specific
            Sess(s).AR(1:SPM.PPM.AR_P,:) = block.ap_mean;
        end

        if SPM.PPM.update_F
            switch SPM.PPM.priors.A
                case 'Robust',
                    Fn=block.F;
                    SPM.PPM.Sess(s).block(z).F=sum(Fn);
                otherwise
                    SPM.PPM.Sess(s).block(z).F = block.F;
                    % Contribution map sums over sessions
                    Fn = spm_vb_Fn(R0Y,block);
            end
            LogEv = LogEv+Fn;
        end

        %-Update regression coefficients
        %---------------------------------------------------------------
        ncols=length(SPM.Sess(s).col);
        beta(SPM.Sess(s).col,:) = block.wk_mean(1:ncols,:);
        if ncols==0
            % Design matrix empty except for constant
            mean_col_index=s;
        else
            mean_col_index=SPM.Sess(nsess).col(end)+s;
        end
        beta(mean_col_index,:) = block.wk_mean(ncols+1,:); % Session mean

        %-Report session-specific noise variances
        %---------------------------------------------------------------
        Sess(s).Hp(1,:)        = sqrt(1./block.mean_lambda');

        %-Store regression coefficient posterior standard deviations
        %---------------------------------------------------------------
        Psd (SPM.Sess(s).col,:) = block.w_dev(1:ncols,:);
        Psd (mean_col_index,:) = block.w_dev(ncols+1,:);

        %-Update contrast variance
        %---------------------------------------------------------------
        if ncon > 0
            for ic=1:ncon,
                CC=SPM.xCon(ic).c;
                % Get relevant columns of contrast
                CC=[CC(SPM.Sess(s).col) ; 0];
                for i=1:nVox,
                    con_var(ic,i)=con_var(ic,i)+CC'*block.w_cov{i}*CC;
                end
            end
        end
        
        switch SPM.PPM.priors.A,
            case 'Robust',
                % Save voxel data where robust model is favoured
                outlier_voxels=find(Fn>0);
                N_outliers=length(outlier_voxels);
                Y_out=R0Y(:,outlier_voxels);
                gamma_out=block.gamma(:,outlier_voxels);
                analysed_xyz=xyz;
                outlier_xyz=analysed_xyz(:,outlier_voxels);

                SPM.PPM.Sess(s).block(z).outlier_voxels=outlier_voxels;
                SPM.PPM.Sess(s).block(z).N_outliers=N_outliers;
                SPM.PPM.Sess(s).block(z).Y_out=Y_out;
                SPM.PPM.Sess(s).block(z).gamma_out=gamma_out;
                SPM.PPM.Sess(s).block(z).outlier_xyz=outlier_xyz;

                block = spm_vb_taylor_R(R0Y,block);
                SPM.PPM.Sess(s).block(z).mean=block.mean;
                SPM.PPM.Sess(s).block(z).N=block.N;
                
                % Prior precision 
                SPM.PPM.Sess(s).block(z).mean_alpha=block.mean_alpha;
                
            otherwise
                % Prior precision
                SPM.PPM.Sess(s).block(z).mean_alpha=block.mean_alpha;
                
                %-Get block-wise Taylor approximation to posterior correlation
                %-------------------------------------------------------
                block = spm_vb_taylor_R(R0Y,block);
                SPM.PPM.Sess(s).block(z).mean=block.mean;
                SPM.PPM.Sess(s).block(z).elapsed_seconds=block.elapsed_seconds;

                %-Save Coefficient RESELS and number of voxels
                %-------------------------------------------------------
                SPM.PPM.Sess(s).block(z).gamma_tot=block.gamma_tot;
                SPM.PPM.Sess(s).block(z).N=block.N;
        end

        %-Save typical structure-specific AR coeffs
        %---------------------------------------------------------------
        if strcmp(SPM.PPM.priors.A,'Discrete')
            SPM.PPM.Sess(s).block(z).as_mean=block.as;
            SPM.PPM.Sess(s).block(z).as_dev=sqrt(1./block.mean_beta);
        end

        clear block;
    end % loop over sessions

    %-Get contrasts
    %-------------------------------------------------------------------
    if ncon > 0
        for ic=1:ncon
            CC=SPM.xCon(ic).c;
            con(ic,:)=CC'*beta;
        end
    end

    %-Append new inmask voxel locations and volumes
    %-------------------------------------------------------------------
    XYZ(:,S + [1:nVox]) = xyz;          %-InMask XYZ voxel coords
    labels(1,S + [1:nVox]) = ones(1,nVox)*z;   %-InMask labels
    S                   = S + nVox;     %-Volume analysed (voxels)

    %-Write conditional beta images
    %-------------------------------------------------------------------
    if z == 1, j  = NaN(xdim,ydim,zdim); end
    for i = 1:nBeta
        if z > 1, j = spm_read_vols(spm_vol(Vbeta(i))); end
        j(Q) = beta(i,:);
        Vbeta(i)  = spm_write_plane(Vbeta(i),j,':'); % this is slow.
        % faster to find and save relevant slices instead of whole volume
    end


    %-Write SD error images
    %-------------------------------------------------------------------
    for s=1:nsess
        if z == 1, j  = NaN(xdim,ydim,zdim); end
        if z > 1, j = spm_read_vols(spm_vol(SPM.PPM.Sess(s).VHp)); end
        j(Q) = Sess(s).Hp(1,:);
        SPM.PPM.Sess(s).VHp = spm_write_plane(SPM.PPM.Sess(s).VHp,j,':');
    end

    %-Write posterior standard-deviation of beta images
    %-------------------------------------------------------------------
    if z == 1, j  = NaN(xdim,ydim,zdim); end
    for i = 1:nPsd
        if z > 1, j = spm_read_vols(spm_vol(VPsd(i))); end
        j(Q) = Psd(i,:);
        VPsd(i) = spm_write_plane(VPsd(i),j,':');
    end

    %-Write AR images
    %-------------------------------------------------------------------
    for s = 1:nsess
        if z == 1, j  = NaN(xdim,ydim,zdim); end
        for i = 1:SPM.PPM.AR_P
            if z > 1, j = spm_read_vols(spm_vol(SPM.PPM.Sess(s).VAR(i))); end
            j(Q) = Sess(s).AR(i,:);
            SPM.PPM.Sess(s).VAR(i) = spm_write_plane(SPM.PPM.Sess(s).VAR(i),j,':');
        end
    end

    %-Write contribution image
    %-------------------------------------------------------------------
    if SPM.PPM.update_F
        if z == 1, j  = NaN(xdim,ydim,zdim); end
        if z > 1, j = spm_read_vols(spm_vol(SPM.PPM.LogEv)); end
        j(Q) = LogEv;
        SPM.PPM.LogEv    = spm_write_plane(SPM.PPM.LogEv,j,':');
    end

    %-Write contrast and contrast SD images
    %-------------------------------------------------------------------
    if ncon > 0
        if z == 1, j  = NaN(xdim,ydim,zdim); end
        for ic=1:ncon
            if z > 1, j = spm_read_vols(spm_vol(SPM.xCon(ic).Vcon)); end
            j(Q) = con(ic,:);
            SPM.xCon(ic).Vcon  = spm_write_plane(SPM.xCon(ic).Vcon,j,':');
        end
        if z == 1, j  = NaN(xdim,ydim,zdim); end
        for ic=1:ncon
            if z > 1, j = spm_read_vols(spm_vol(SPM.PPM.Vcon_sd(ic))); end
            j(Q)  = sqrt(con_var(ic,:));
            SPM.PPM.Vcon_sd(ic) = spm_write_plane(SPM.PPM.Vcon_sd(ic),j,':');
        end
    end

    if SPM.PPM.window
        %-Report progress
        %---------------------------------------------------------------
        spm_progress_bar('Set',100*z/nLb);
    end

end % (for z = 1:nLb)

%-Done!
%-----------------------------------------------------------------------
fprintf('\n')                                                           %-#
if SPM.PPM.window
    spm_progress_bar('Clear')
    spm('Pointer','Arrow');
end

%=======================================================================
% - P O S T   E S T I M A T I O N
%=======================================================================

if S == 0, warning('No inmask voxels - empty analysis!'), end

%-Create 1st contrast for 'effects of interest' (all if not specified)
%=======================================================================
Fcname     = 'effects of interest';
try
    iX0    = [xX.iB xX.iG];
catch
    iX0    = [];
end

xX.xKXs = spm_sp('Set',spm_filter(xX.K,xX.X));      % ** Not Whitened **
xX.erdf = size(xX.X,1); % Just set to number of scans so, when 
                        % we assess the results, spm_getSPM is happy
xX.W= eye(size(xX.X,1)); % Set whitening matrix to identity -
                         % we must set it to keep spm_graph happy
                         
xCon       = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);

%-Compute scaled design matrix for display purposes
%-----------------------------------------------------------------------
xX.nKX     = spm_DesMtx('sca',xX.xKXs.X,xX.name);

%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                    %-#

%-place fields in SPM
%-----------------------------------------------------------------------
SPM.xVol.XYZ   = XYZ(:,1:S);    %-InMask XYZ coords (voxels)
SPM.xVol.labels= labels(:,1:S); %-InMask labels (voxels)
SPM.xVol.M     = M;             %-voxels -> mm
SPM.xVol.iM    = inv(M);        %-mm -> voxels
SPM.xVol.DIM   = DIM;           %-image dimensions
SPM.xVol.S     = S;             %-Volume (voxels)
SPM.xVol.R     = 100;           % Set R - number of RESELS - to arbitrary value
                                % as, if R not set, SPM will think model has not 
                                % been estimated
SPM.xVol.FWHM  = 10;            % Set to arbitrary value so spm_getSPM is happy
                              
SPM.VCbeta     = Vbeta;         % Filenames - parameters
SPM.VPsd       = VPsd;          % Filenames - hyperparameters
SPM.VM         = VM;            %-Filehandle - Mask

SPM.PPM.Gamma  = 1;             % Default threshold for effect size (1 per cent)

SPM.xX         = xX;            %-design structure
SPM.xM         = xM;            %-mask structure

% Copy contrast structure 
SPM.PPM.xCon   = SPM.xCon;
for i=1:length(SPM.PPM.xCon)
    SPM.PPM.xCon(i).PSTAT='T';
end

% Add pointer to RPV image file so that spm_list works
SPM.xVol.VRpv=[];

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat', 'SPM', fmt);

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#

if SPM.PPM.window
    %===================================================================
    %- E N D: Cleanup GUI
    %===================================================================
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    fprintf('%-40s: %30s\n','Completed',spm('time'))                    %-#
    fprintf('...use the results section for assessment\n\n')            %-#
end
