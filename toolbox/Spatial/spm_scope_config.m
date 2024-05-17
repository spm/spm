function scope = spm_scope_config
% SPM Configuration file for SCOPE distortion correction
%__________________________________________________________________________

% Nicole Labra Avila
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% SCOPE
%--------------------------------------------------------------------------
scope      = cfg_exbranch;
scope.tag  = 'scope';
scope.name = 'SCOPE';
scope.val  = @scope_cfg;
scope.prog = @(job)spm_run_scope('run',job);
scope.vout = @(job)spm_run_scope('vout',job);
scope.help = {[...
'Utility to correct susceptibility distortions in EPI images using ',...
'a re-implementation of FSL''s topup (J.L.R. Andersson, S. Skare, J. Ashburner, 2003. ',...
'How to correct susceptibility distortions in spin-echo echo-planar ',...
'images: application to diffusion tensor imaging). This implementation ',...
'requires two EPI images acquired with opposite traversals along the phase-encode direction, ',...
'which is assumed to be in the y-direction (in voxel-space).']};

%==========================================================================
function varargout = scope_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%==========================================================================
%--------------------------------------------------------------------------
% Same blip direction volume(s)
%--------------------------------------------------------------------------
vol_same         = cfg_files;
vol_same.tag     = 'vol1';
vol_same.name    = 'Same PE direction image(s)';
vol_same.filter  = 'image';
vol_same.ufilter = '.*';
vol_same.num     = [1 inf];
vol_same.help    = {...
['Select an image, or set of images, where the acquisition traversed the phase-encode '...,
 'direction of K-space in the *same* way as the fMRI data to be corrected ',...
 '(and could even be the same data as the fMRI to be corrected). ',...
 'Note that if multiple images are specified, then they will be motion corrected and ',...
 'their average used for registering the phase-encoding direction reversed images. ',...
 'If these images were acquired before the oposite PE images, then it is advised to select the ',...
 'final volume before the others, otherwise simply select them in the usual chronological order.']};
vol_same.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Opposite blip direction volume(s)
%--------------------------------------------------------------------------
vol_oppo         = cfg_files;
vol_oppo.tag     = 'vol2';
vol_oppo.name    = 'Opposite PE direction image(s)';
vol_oppo.filter  = 'image';
vol_oppo.ufilter = '.*';
vol_oppo.num     = [1 inf];
vol_oppo.help    = {...
['Select an image, or set of images, where the acquisition traversed the phase-encode ',...
 'direction of K-space in the *opposite* way to the fMRI data to be corrected. ',...
 'Note that if multiple images are specified, then they will be motion corrected and ',...
 'their average used for registering the phase-encoding direction reversed images .',...
 'If these images were acquired before the same PE images, then it is advised to select the ',...
 'final volume before the others, otherwise simply select them in the usual chronological order.']}; 
vol_oppo.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% fwhm values
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.val     = {[8 4 2 1 0]};
fwhm.strtype = 'e';
fwhm.num     = [1 Inf];
fwhm.help    = {[...
'Spatial scales (expressed as FWHM of Gaussian kernels used to smooth ',...
'input images during SCOPE estimation.']};

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
rinterp.name    = 'Interpolation';
rinterp.help    = {[
'Degree of B-spline (from 0 to 7) along different dimensions ' ...
'(see ``spm_diffeo``).']};
rinterp.labels = {
                  'Linear              '
                  '2nd Degree B-spline '
                  '3rd Degree B-Spline'
                  '4th Degree B-Spline'
                  '5th Degree B-Spline '
                  '6th Degree B-Spline'
                  '7th Degree B-Spline'
}';
rinterp.values = {1,2,3,4,5,6,7};
rinterp.val    = {2};

%--------------------------------------------------------------------------
% prefix VDM Filename Prefix
% Option for Jacobian scaling
%--------------------------------------------------------------------------
jac         = cfg_menu;
jac.tag     = 'jac';
jac.name    = 'Jacobian scaling';
jac.val     = {1};
jac.help    = {
    'Option to include Jacobian scaling in the registration model.'
    ['It includes in the process the changes of intensities due to' ...
    ' stretching and compression.']
    }';
jac.labels  = {
              'Yes'
              'No'
}';
jac.values  = {1 0};

%--------------------------------------------------------------------------
% prefix VDM Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'vdm Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the voxel displacement map (vdm) file. Default prefix is ``vdm5_``.'};
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
'All created files (voxel displacement map and unwarped images) are written ',...
'in the specified directory. The voxel displacement map is saved to disk as ',...
'a vdm file (``vdm5_*.nii``)']};

[cfg,varargout{1}] = deal({vol_same,vol_oppo,fwhm,reg,rinterp,jac,prefix,outdir});

%==========================================================================
function out = spm_run_scope(cmd, job)

switch lower(cmd)
    case 'run'
        vdm               = spm_scope(job.vol1,job.vol2,job.fwhm,job.reg, ...
                                      job.rinterp,job.jac,job.prefix,job.outdir{1});
        out.vdmfile       = {vdm.dat.fname};

    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Voxel displacement map';
        out(1).src_output = substruct('.','vdmfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end


