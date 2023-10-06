function cfg = spm_TVdenoise_config
% SPM Configuration file for total variation denoising
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% topup Topup
%--------------------------------------------------------------------------
cfg      = cfg_exbranch;
cfg.tag  = 'denoise';
cfg.name = 'Total-variation denoising';
cfg.val  = @denoise_cfg;
cfg.prog = @(job)spm_run_denoise('run',job);
cfg.vout = @(job)spm_run_denoise('vout',job);
cfg.help = {'Total-variation (multi-channel) denoising of brain MRI.'};

%==========================================================================
function varargout = denoise_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%--------------------------------------------------------------------------
% vols Volumes
%--------------------------------------------------------------------------
vols         = cfg_files;
vols.tag     = 'data';
vols.name    = 'Volumes';
vols.help    = {
'Select scans from this channel for processing. If multiple channels are used (eg PDw & T2w), then the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..).'
               }';
vols.filter  = 'image';
vols.ufilter = '.*';
vols.num     = [1 Inf];
vols.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% data Data
%--------------------------------------------------------------------------
data         = cfg_repeat;
data.tag     = 'data';
data.name    = 'Data';
data.val     = {vols};
data.help    = {
'Specify the number of different channels (for multi-spectral classification). If you have scans of different contrasts for each of the subjects, then it is possible to combine the information from them in order to improve the segmentation accuracy. Note that only the first channel of data is used for the initial affine registration with the tissue probability maps.'
               }';
data.values  = {vols};
data.num     = [1 Inf];

lambda        = cfg_menu;
lambda.tag    = 'lambda';
lambda.name   = 'Denoising strength';
lambda.values = {10, 30, 100, 300, 1000};
lambda.labels = {'  10: Weak', '  30: Reasonable', ' 100: Strong',' 300: Very strong', '1000: For testing'};
lambda.val    = {30};
lambda.help   = {'Vary the strength of denoising.'};

nit        = cfg_menu;
nit.tag    = 'nit';
nit.name   = 'Iterations';
nit.values = {10, 30, 100, 300, 1000};
nit.labels = {'  10: Fastest/poorest', '  30: Fast/poor', ' 100: reasonable',' 300: Slow/good', '1000: For testing'};
nit.val    = {100};
nit.help   = {'Number of denoising relaxation iterations.'};

dev        = cfg_menu;
dev.tag    = 'dev';
dev.name   = 'Device';
dev.values = {0, 1};
dev.labels = {'CPU', 'GPU'};
dev.val    = {0};
dev.help   = {'Run on CPU/GPU.'};

[cfg,varargout{1}] = deal({data lambda nit dev});

