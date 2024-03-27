function scope = spm_scope_config
% SPM Configuration file for Blip updown distortion correction
%__________________________________________________________________________

% Nicole Labra Avila
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% Blip Updown
%--------------------------------------------------------------------------
scope      = cfg_exbranch;
scope.tag  = 'scope';
scope.name = 'SCOPE';
scope.val  = @scope_cfg;
scope.help = {'Correct susceptibility distortions using Blip Updown.'};
scope.prog = @(job)spm_run_scope('run',job);
scope.vout = @(job)spm_run_scope('vout',job);
scope.help = {[...
'Utility to correct susceptibility distortions in EPI images using ',...
'a re-implementation of FSL''s topup (J.L.R. Andersson, S. Skare, J. Ashburner, 2003. ',...
'How to correct susceptibility distortions in spin-echo echo-planar ',...
'images: application to diffusion tensor imaging). This implementation ',...
'requires two EPI images acquired with opposed phase-encode blips in ',...
'the y-direction (for the moment).']};

%==========================================================================
function varargout = scope_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%==========================================================================
%--------------------------------------------------------------------------
% Blip up volume(s)
%--------------------------------------------------------------------------
volbup         = cfg_files;
volbup.tag     = 'volbup';
volbup.name    = 'Blip reversed data with positive polarity (.nii)';
volbup.filter  = 'image';
volbup.ufilter = '.*';
volbup.num     = [1 inf];
volbup.help    = {'Select a blip reversed image with positive phase-encoding polarity.'};
volbup.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Blip down volume(s)
%--------------------------------------------------------------------------
volbdown         = cfg_files;
volbdown.tag     = 'volbdown';
volbdown.name    = 'Blip reversed data with negative polarity (.nii)';
volbdown.filter  = 'image';
volbdown.ufilter = '.*';
volbdown.num     = [1 inf];
volbdown.help    = {'Select an blip reversed image with negative phase-encoding polarity.'};
volbdown.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Data Volumes
%--------------------------------------------------------------------------
data         = cfg_branch;
data.tag     = 'data';
data.name    = 'Volumes';
data.val     = {volbup volbdown};
data.help    = {[....*
'Selectthe images with opposite phase-encoding polarities. The first one(s) ',...
'must have positive polarity and the second one(s) negative polarity.']};

%--------------------------------------------------------------------------
% Acquisition of images order
%--------------------------------------------------------------------------
acqorder         = cfg_menu;
acqorder.tag     = 'acqorder';
acqorder.name    = 'Acquisition order of images';
acqorder.help    = {[
'Order of acquisition of blip-reversed images ']};
acqorder.labels = {
                  'Positive polarity first'
                  'Negative polarity first'
}';
acqorder.values = {0,1};
acqorder.val    = {0};

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
'input images during Blip Updown estimation.']};

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

[cfg,varargout{1}] = deal({data,acqorder,fwhm,reg,rinterp,jac,prefix,outdir});


%==========================================================================
function out = spm_run_scope(cmd, job)

switch lower(cmd)
    case 'run'
        vdm               = spm_scope(job.data,job.acqorder,job.fwhm,job.reg, ...
                            job.rinterp,job.jac,job.prefix,job.outdir{1});
        out.vdmfile       = {vdm.dat.fname};

    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Voxel displacement map';
        out(1).src_output = substruct('.','vdmfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

end


