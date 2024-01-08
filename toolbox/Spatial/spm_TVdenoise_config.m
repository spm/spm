function cfg = spm_TVdenoise_config
% SPM Configuration file for total variation denoising
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% denoise Total-variation denoising
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
'Specify the number of different image channels. If you have scans of different contrasts for each of the subjects, then it is possible to combine the information from them in a way that might improve the denoising.'
               }';
data.values  = {vols};
data.num     = [1 Inf];

lambda         = cfg_entry;
lambda.tag     = 'lambda';
lambda.name    = 'Denoising strength';
lambda.strtype = 'r';
lambda.num     = [1 1];
lambda.val     = {0.03};
lambda.help    = {'Vary the strength of denoising. This setting may need tweaking to obtain the best results. Note that the standard deviation over the field of view is used to adjust this setting.'};

nit        = cfg_menu;
nit.tag    = 'nit';
nit.name   = 'Iterations';
nit.values = {10, 30, 100, 300, 1000};
nit.labels = {'  10: Fastest/poorest', '  30: Fast/poor', ' 100: reasonable',' 300: Slow/good', '1000: For testing'};
nit.val    = {100};
nit.help   = {'Number of denoising relaxation iterations.'};

% Better disable this option because SPM is supposed to not require
% additional toolbox licenses.
% I'm not sure what the proper way to test if a license for the Distributed
% Computing Toolbox exists.
if false % license('test','distrib_computing_toolbox')
    device        = cfg_menu;
    device.tag    = 'device';
    device.name   = 'Device';
    device.values = {'cpu','gpu'};
    device.labels = {'CPU', 'GPU'};
    device.val    = {'cpu'};
    device.help   = {'Run on CPU/GPU. Note that using gpuArrays requires MATLAB''s Distributed Computing Toolbox.'};
else
    device          = cfg_const;
    device.tag      = 'device';
    device.val      = {'cpu'};
    device.hidden   = true;
end

[cfg,varargout{1}] = deal({data lambda nit device});

