function topup = spm_topup_config
% SPM Configuration file for Topup
%__________________________________________________________________________

% Nicole Labra Avila
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% topup Topup
%--------------------------------------------------------------------------
topup      = cfg_exbranch;
topup.tag  = 'topup';
topup.name = 'Topup';
topup.val  = @topup_cfg;
topup.help = {'Correct susceptibility distortions using topup.'};
topup.prog = @(job)spm_run_topup('run',job);
topup.vout = @(job)spm_run_topup('vout',job);
topup.help = {[...
'Utility to correct susceptibility distortions in EPI images using ',...
'the topup algorithm (J.L.R. Andersson, S. Skare, J. Ashburner, 2003. ',...
'How to correct susceptibility distortions in spin-echo echo-planar ',... 
'images: application to diffusion tensor imaging). This implementation ',...
'requires two EPI images acquired with opposed phase-encode blips in ',...
'the y-direction (for the moment).']};

%==========================================================================
function varargout = topup_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%==========================================================================
%--------------------------------------------------------------------------
% Blip up volume(s)
%--------------------------------------------------------------------------
volbup         = cfg_files;
volbup.tag     = 'volbup';
volbup.name    = 'Blip up volume (.nii)';
volbup.filter  = 'image';
volbup.ufilter = '.*';
volbup.num     = [0 1];
volbup.val     = {''};
volbup.help    = {'Select an image with a blip up phase-encoding polarity.'};
volbup.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Blip down volume(s)
%--------------------------------------------------------------------------
volbdown         = cfg_files;
volbdown.tag     = 'volbdown';
volbdown.name    = 'Blip down volume (.nii)';
volbdown.filter  = 'image';
volbdown.ufilter = '.*';
volbdown.num     = [0 1];
volbdown.val     = {''};
volbdown.help    = {'Select an image with a blip down phase-encoding polarity.'};
volbdown.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% data Volumes
%--------------------------------------------------------------------------
data         = cfg_branch;
data.tag     = 'data';
data.name    = 'Volumes';
data.val     = {volbup volbdown};
data.help    = {[....*
'Select two images with opposite phase-encoding polarities. The first one ',...
'a blip up and the second one a blip down.']};

%--------------------------------------------------------------------------
% fwhm values
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.val     = {[8 4 2 1 0.1]};
fwhm.strtype = 'e';
fwhm.num     = [1 Inf];
fwhm.help    = {[...
'Spatial scales (expressed as FWHM of Gaussian kernels used to smooth ',...
'input images during topup estimation.']};

%--------------------------------------------------------------------------
% Regularisation .*
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Regularisation';
reg.val     = {[0 10 100]}; % Default values 
reg.strtype = 'e';
reg.num     = [1 3];
reg.help    = {...
'Regularisation settings (see spm_field).',...
'The three expected values are:',...
'    1. Penalty on absolute values.',...
'    2. Penalty on the *membrane energy*.',...
'    3. Penalty on the *bending energy*.'};

%--------------------------------------------------------------------------
% Degree of B-spline  
%--------------------------------------------------------------------------
rinterp         = cfg_menu;
rinterp.tag     = 'rinterp';
rinterp.name        = 'Interpolation';
rinterp.val     = {[1 1 1]};         % Default values 
rinterp.help    = {[
'Degree of B-spline (from 0 to 7) along different dimensions ' ...
'(see ``spm_diffeo``).']};
rinterp.labels = {
                  '2nd Degree B-spline '
                  '3rd Degree B-Spline'
                  '4th Degree B-Spline'
                  '5th Degree B-Spline '
                  '6th Degree B-Spline'
                  '7th Degree B-Spline'
}';
rinterp.values = {[0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1] [1 1 1]};

%--------------------------------------------------------------------------
% prefix VDM Filename Prefix
% Option for refine topup
%--------------------------------------------------------------------------
rt         = cfg_menu;
rt.tag     = 'rt';
rt.name    = 'Refine topup';
rt.val     = {1};
rt.help    = {
    'Option to include refine topup after estimating the fields.'
    ['It includes in the process the changes of intensities due to' ...
    ' stretching and compression.']
    }';
rt.labels  = {
              'Yes'
              'No'
}';
rt.values  = {1 0};

%--------------------------------------------------------------------------
% prefix VDM Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'VDM Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the VDM files. Default prefix is ``vdm5_``.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'vdm5'};


%--------------------------------------------------------------------------
% Save VDM 
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.name    = 'Output directory';
outdir.tag     = 'outdir';
outdir.val{1}  = {''};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
outdir.help    = {[...
'All created files (deformation fields and unwarped images) are written ',...
'in the specified directory. The deformation field is saved to disk as ',...
'a vdm file (``vdm5_*.nii``)']};

[cfg,varargout{1}] = deal({data,fwhm,reg,rinterp,rt,prefix,outdir});


%==========================================================================
function out = spm_run_topup(cmd, job)

switch lower(cmd)
    case 'run'
        VDM               = spm_topup(job.data.volbup{1},job.data.volbdown{1}, ...
                            job.fwhm,job.reg,job.rinterp,job.rt,job.prefix, ...
                            job.outdir{1});
        out.vdmfile       = {VDM.dat.fname};
        
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Voxel displacement map';
        out(1).src_output = substruct('.','vdmfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

end
