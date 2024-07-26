function [SPM] = spm_spm_Bayes(SPM)
% Conditional parameter estimation of a General Linear Model
% FORMAT [SPM] = spm_spm_Bayes(SPM)
%__________________________________________________________________________
%
% spm_spm_Bayes returns to voxels identified by spm_spm (ML parameter
% estimation) to get conditional parameter estimates and ReML hyper-
% parameter estimates.  These estimates use prior covariances, on the
% parameters, from empirical Bayes.  These PEB prior variances come from
% the hierarchical model that obtains by considering voxels as providing a
% second level.  Put simply, the variance in parameters, over voxels, is
% used as a prior variance from the point of view of any one voxel. The
% error covariance hyperparameters are re-estimated in the light of these
% priors.  The approach adopted is essentially a fully Bayesian analysis at
% each voxel, using empirical Bayesian prior variance estimators over
% voxels.
%
% Each separable partition (i.e. session) is assigned its own
% hyperparameter but within session covariance components are lumped
% together, using their relative expectations over voxels.  This makes
% things much more computationally efficient and avoids inefficient
% voxel-specific multiple hyperparameter estimates.
%
% spm_spm_Bayes adds the following fields to SPM:
%
%                           ----------------
%
%
%   SPM.PPM.l      = session-specific hyperparameter means
%   SPM.PPM.Cb     = empirical prior parameter covariances
%   SPM.PPM.C      = conditional covariances of parameters
%   SPM.PPM.dC{i}  = dC/dl;
%   SPM.PPM.ddC{i} = ddC/dldl
%
% The derivatives are used to compute the conditional variance of various
% contrasts in spm_getSPM, using a Taylor expansion about the hyperparameter
% means.
%
%
%                           ----------------
%
%   SPM.VCbeta     - Handles of conditional parameter estimates
%   SPM.VHp        - Handles of hyperparameter estimates
%
%                           ----------------
%
% Cbeta_????.<ext>                     - conditional parameter images
% These are 32-bit (float) images of the conditional estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.<ext>) are given
% value NaN.
%
%                           ----------------
%
% CHp_????.<ext>              - error covariance hyperparameter images
% This is a 64-bit (double) image of the ReML error variance estimate.
% for each separable partition (Session).  Voxels outside the analysis 
% mask are given value NaN.
%__________________________________________________________________________
% 
% For single subject fMRI analysis there is an alternative function
% using voxel-wise GLM-AR models that are spatially regularised
% using the VB framework. This is implemented using spm_spm_vb.m.
%__________________________________________________________________________

% Karl Friston, Peter Zeidman
% Copyright (C) 2002-2022 Wellcome Centre for Human Neuroimaging

% Whether to use parallel computing (this will be moved to spm_defaults in
% a future relase)
use_parfor = false;

%-Say hello
Finter = spm('FigName','Stats: Bayesian estimation...');

% Select SPM.mat & change directory
%--------------------------------------------------------------------------
if ~nargin
    [Pf, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd = spm_file(Pf,'fpath');
    load(fullfile(swd,'SPM.mat'))
    cd(swd)
end

% Unpack inputs
%--------------------------------------------------------------------------
try
    M    = SPM.xVol.M;
    DIM  = SPM.xVol.DIM;
    xdim = DIM(1); ydim = DIM(2); zdim = DIM(3);
    XYZ  = SPM.xVol.XYZ;
catch
    helpdlg({   'Please do a ML estimation first.',...
            'This identifies the voxels to analyse.'});
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
end

%-Initialise output images
%--------------------------------------------------------------------------
fprintf('%-40s: %30s','Output images','...initialising')

%-Initialise conditional estimate image files
xX             = SPM.xX;
[nScan,nBeta]  = size(xX.X);
Vbeta(1:nBeta) = deal(struct(...
            'fname',   [],...
            'dim',     DIM',...
            'dt',      [spm_type('float32'), spm_platform('bigend')],...
            'mat',     M,...
            'pinfo',   [1 0 0]',...
            'descrip', ''));
for i = 1:nBeta
    Vbeta(i).fname   = [sprintf('Cbeta_%04d',i) spm_file_ext];
    Vbeta(i).descrip = sprintf('Cond. beta (%04d) - %s',i,xX.name{i});
    spm_unlink(Vbeta(i).fname)
end
Vbeta = spm_create_vol(Vbeta);

%-Initialise ReML hyperparameter image files
try
    nHp       = length(SPM.nscan);
catch
    nHp       = nScan;
    SPM.nscan = nScan;
end

% Number of separable partitions
s = nHp;

VHp(1:nHp)    = deal(struct(...
            'fname',   [],...
            'dim',     DIM',...
            'dt',      [spm_type('float64'), spm_platform('bigend')],...
            'mat',     M,...
            'pinfo',   [1 0 0]',...
            'descrip', ''));
for i = 1:nHp
    VHp(i).fname   = [sprintf('Hp_%04d',i) spm_file_ext];
    VHp(i).descrip = sprintf('Hyperparameter (%04d)',i);
    spm_unlink(VHp(i).fname)
end
VHp   = spm_create_vol(VHp);

% Estimate covariances using whole-brain data
% -------------------------------------------------------------------------
fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimating data covariance')
SPM.xVi.CY = spm_spm_Bayes_CY(SPM);

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...estimating priors')
sP = spm_spm_Bayes_specify(SPM);

% Fit model and write images
% -------------------------------------------------------------------------

% Report
spm_progress_bar('Init',100,'Bayesian estimation','');
spm('Pointer','Watch')

if use_parfor
    % Create or get parallel pool
    pool = gcp;
    
    % Get the number of workers
    nw = pool.NumWorkers;
else
    nw  = 1;
end

% Get maximum amount of data processed at a time (bytes)
MAXMEM = spm_get_defaults('stats.maxmem');

% Max block size - the number of voxel-wise vectors that will fit in
% memory per worker (based on 8 byte, i.e. 64 bit, images)
max_blksz = ceil(MAXMEM/8/nScan/nw);

% Sum of hyperparameters
SHp = 0;

VY = SPM.xY.VY;

% Loop over planes
for z = 1:zdim

    U       = find(XYZ(3,:) == z);    
    CrBl    = zeros(nBeta,length(U)); % conditional parameter estimates
    CrHp    = zeros(nHp,  length(U)); % ReML hyperparameter estimates

    % Set maximum size of one block in voxels
    if use_parfor && ~isempty(U)
        % Equal distribution of voxels among workers
        blksz = ceil(length(U) / nw);
        blksz = min(blksz, max_blksz); 
    else
        % All voxels in one block (if possible)
        blksz = max_blksz;
    end
    
    % Number of blocks
    nbch = ceil(length(U)/blksz);    

    % Construct list of voxel indices in U and coordinates for each block    
    idx = cell(1, nbch);
    xyz = cell(1, nbch);
    for bch = 1:nbch
        I     = (1:blksz) + (bch - 1)*blksz;
        I     = I(I <= length(U));

        idx{bch} = I;
        xyz{bch} = XYZ(:,U(I));
    end
    
    beta_blocks = cell(1,nbch);
    Hp_blocks   = cell(1,nbch);
    
    % Loop over blocks (bunches of lines / planks)
    if use_parfor
        parfor bch = 1:nbch
            % Get response variable
            Y = spm_get_data(VY,xyz{bch});
            
            % Invert
            [beta,Hp] = invert_glm(Y,sP);
            
            % Store
            beta_blocks{bch} = beta;
            Hp_blocks{bch} = Hp;
        end
    else
        for bch = 1:nbch
            % Get response variable
            Y = spm_get_data(VY,xyz{bch});
            
            % Invert
            [beta,Hp] = invert_glm(Y,sP);
            
            % Store
            beta_blocks{bch} = beta;
            Hp_blocks{bch} = Hp;
        end
    end
    
    if nbch > 0
        % Concatenate blocks
        idx = cell2mat(idx);
        beta_blocks = cell2mat(beta_blocks);
        Hp_blocks   = cell2mat(Hp_blocks);

        % Store in image plane
        CrBl(:,idx) = beta_blocks;
        CrHp(:,idx) = Hp_blocks;

        % Sum hyperparameters    
        SHp = SHp + sum(Hp_blocks,2);
    end

    % Write plane to conditional beta images
    for i = 1:nBeta
        tmp       = sparse(XYZ(1,U),XYZ(2,U),CrBl(i,:),xdim,ydim);
        tmp(~tmp) = NaN;
        Vbeta(i)  = spm_write_plane(Vbeta(i),tmp,z);
    end

    % Write plane to hyperparameter images
    for i = 1:nHp
        tmp       = sparse(XYZ(1,U),XYZ(2,U),CrHp(i,:),xdim,ydim);
        tmp(~tmp) = NaN;
        VHp(i)    = spm_write_plane(VHp(i),tmp,z);
    end

    % Report progress
    spm_progress_bar('Set',100*(z - 1)/zdim);

end % (for z = 1:zdim)

fprintf('\n')
spm_progress_bar('Clear')

%==========================================================================
% - P O S T   E S T I M A T I O N
%==========================================================================

% Taylor expansion for conditional covariance
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Non-sphericity','...REML estimation')          %-#

% expansion point (mean hyperparameters)
l     = SHp/SPM.xVol.S;

% change in conditional covariance w.r.t. hyperparameters
%--------------------------------------------------------------------------
n     = size(xX.X,2);
PPM.l = l;
for i = 1:s
    PPM.dC{i}  = sparse(n,n);
    PPM.ddC{i} = sparse(n,n);
end
for i = 1:s

    P     = sP(i).P;
    u     = sP(i).u;
    v     = sP(i).v;

    % derivatives of conditional covariance w.r.t. hyperparameters
    %----------------------------------------------------------------------
    d     = P{1}.X'*inv(P{1}.C{1})*P{1}.X;
    Cby   = inv(d/l(i) + inv(P{2}.C));
    d     = d*Cby;
    dC    = Cby*d/(l(i)^2);
    ddC   = 2*(dC/(l(i)^2) - Cby/(l(i)^3))*d;

    % place in output structure
    %----------------------------------------------------------------------
    j               = 1:length(v);
    PPM.Cb(v,v)     = P{2}.C(j,j);
    PPM.Cby(v,v)    = Cby(j,j);
    PPM.dC{i}(v,v)  = dC(j,j);
    PPM.ddC{i}(v,v) = ddC(j,j);        
    
end

%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing')

%-Save analysis parameters in SPM.mat file
%--------------------------------------------------------------------------
SPM.VCbeta = Vbeta;         % Filenames - parameters
SPM.VHp    = VHp;           % Filenames - hyperparameters
SPM.PPM    = PPM;           % PPM structure

fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat', 'SPM', fmt);

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done') 

% Clean up GUI
%--------------------------------------------------------------------------
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))

% -------------------------------------------------------------------------
function [beta,Hp] = invert_glm(Y,sP)
% Invert GLM(s) for one or more voxels
%
% Y  - voxel-wise data [nVol x nVox]
% sP - model spec structure for spm_peb.m
%
% beta - estimated regression parameters [nBeta x nVox]
% Hp   - estimated hyperparameters [nHp x nVox]

nBeta = size(sP.P{1}.X,2); % Number of parameters
nVox  = size(Y,2);         % Number of voxels
nHp   = length(sP);        % Number of hyperparameters/separable partitions

% Conditional estimates (per partition, per voxel)
beta  = zeros(nBeta,nVox);
Hp    = zeros(nHp,  nVox);
for j = 1:nHp
    P     = sP(j).P; % spm_peb model structure
    u     = sP(j).u; % design matrix rows
    v     = sP(j).v; % design matrix columns
    for i = 1:nVox
        C         = spm_PEB(Y(u,i),P,false,true);
        beta(v,i) = C{2}.E(1:length(v));
        Hp(j,i)   = C{1}.h;
    end
end