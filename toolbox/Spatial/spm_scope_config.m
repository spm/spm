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
 'If these images were acquired before the opposite PE images, then it is advised to select the ',...
 'final volume before the others, otherwise simply select them in the usual chronological order.']};
vol_same.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Same blip direction phase volume(s)
%--------------------------------------------------------------------------
vol_same_pha         = cfg_files;
vol_same_pha.tag     = 'vol1_pha';
vol_same_pha.name    = 'Same PE direction phase image(s)';
vol_same_pha.val{1}  = {};
vol_same_pha.filter  = 'image';
vol_same_pha.ufilter = '.*';
vol_same_pha.num     = [0 inf];
vol_same_pha.help    = {...
['Select a phase image, or set of images, where the acquisition traversed the phase-encode '...,
 'direction of K-space in the *same* way as the fMRI data to be corrected ',...
 'If these images were acquired before the opposite PE images, then it is advised to select the ',...
 'final volume before the others, otherwise simply select them in the usual chronological order.']};

%--------------------------------------------------------------------------
% Total readout time in seconds
%--------------------------------------------------------------------------
t_readout         = cfg_entry;
t_readout.tag     = 't_readout';
t_readout.name    = 'Total readout time';
t_readout.val     = {0};
t_readout.strtype = 'e';
t_readout.num     = [0 1];
t_readout.help    = {...
['Specify EPI total readout time in seconds, which is defined as '...
'total_readout_time = matrix_size_in_PE_direction x nominal_echo_spacing(i.e.2xramp_time+1xflat_top_time of a trapezoid EPI readout) / inplane_acceleration / inplane_segments']};

%--------------------------------------------------------------------------
% Echo time in seconds
%--------------------------------------------------------------------------
echo_time         = cfg_entry;
echo_time.tag     = 'TE';
echo_time.name    = 'Echo time';
echo_time.val     = {0};
echo_time.strtype = 'e';
echo_time.num     = [0 1];
echo_time.help    = {'Specify EPI echo time in seconds'};

%--------------------------------------------------------------------------
% Phase encoding direction
%--------------------------------------------------------------------------
PE_dir         = cfg_menu;
PE_dir.tag     = 'PE_dir';
PE_dir.name    = 'fMRI phase polarity';
PE_dir.val     = {1};
PE_dir.labels  = {
              '-1'
              '1'
                }';
PE_dir.values     = {-1 1};
PE_dir.help = {
    'Depending on the scanner vendor, image reconstruction algorithm, and phase-encoding direction,'
    'the phase image may need to be multiplied by −1 to ensure correct distortion correction.'
    'Using the wrong sign will exacerbate distortions rather than correct them.'
    } ;

%--------------------------------------------------------------------------
% Phase-based initial VDM estimate calculation branch
%--------------------------------------------------------------------------

pha_branch = cfg_branch;
pha_branch.tag = 'pha_branch';
pha_branch.name = 'Calculate initial estimate';
pha_branch.val = {vol_same_pha, t_readout, echo_time, PE_dir};
pha_branch.help =  {'Optional: Provide the phase images and sequence parameters required to compute the initial phase-based VDM estimate.'};
pha_branch.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Voxel Displacement Map prior, i.e. starting estimate
%--------------------------------------------------------------------------
vdm_prior         = cfg_files;
vdm_prior.tag     = 'VDMprior';
vdm_prior.name    = 'Load initial estimate';
vdm_prior.val     = {''};
vdm_prior.filter  = 'image';
vdm_prior.ufilter = '.*';
vdm_prior.num     = [0 1];
vdm_prior.help    = {'Select an image containing the initial estimate of the voxel displacement map. '}; 
vdm_prior.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% No Voxel Displacement Map prior, i.e. no starting estimate
%--------------------------------------------------------------------------
no_vdm_prior         = cfg_const;
no_vdm_prior.tag     = 'no_vdm_prior';
no_vdm_prior.name    = 'No initial estimate';
no_vdm_prior.val     = {[]};

%--------------------------------------------------------------------------
% Voxel Displacement Map prior, i.e. starting estimate
%--------------------------------------------------------------------------
vdm_prior_select    = cfg_choice;
vdm_prior_select.tag    = 'vdm_prior_select';
vdm_prior_select.name   = 'VDM - initial estimate';
vdm_prior_select.values = {no_vdm_prior, pha_branch, vdm_prior} ;
vdm_prior_select.val    ={no_vdm_prior} ;
vdm_prior_select.help = {
    'Optional: Compute or provide initial Voxel Displacement Map (VDM) estimate.'
    'Calculating initial VDM estimate can improve distortion correction accuracy especially in highly-distorted scans.'
    'WARNING: On first use, "Calculate initial estimate" feature will automatically download julia and add MriResearchTools julia source code'
    } ;

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
jac.val     = {0};
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
prefix.val     = {'vdm5_'};


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

[cfg,varargout{1}] = deal({vol_oppo, vol_same, vdm_prior_select, fwhm, reg, rinterp, jac, prefix, outdir});

%==========================================================================

function out = spm_run_scope(cmd, job)

switch lower(cmd)
    case 'run'        
        if isfield(job.vdm_prior_select, 'VDMprior')
            vdm               = spm_scope(job.vol1,job.vol2, job.fwhm,job.reg, ...
                                      job.rinterp,job.jac,job.prefix,job.outdir{1}, job.vdm_prior_select.VDMprior);
        elseif isfield(job.vdm_prior_select.pha_branch, 'vol1_pha')
            vdm               = spm_scope(job.vol1,job.vol2, job.fwhm,job.reg, ...
                                      job.rinterp,job.jac,job.prefix,job.outdir{1}, [],job.vdm_prior_select.pha_branch.vol1_pha, job.vdm_prior_select.pha_branch.t_readout, job.vdm_prior_select.pha_branch.TE, job.vdm_prior_select.pha_branch.PE_dir);
        end
        out.vdmfile       = {vdm.dat.fname};

    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Voxel displacement map';
        out(1).src_output = substruct('.','vdmfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end


