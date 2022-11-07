function topup = spm_topup_config
% SPM Configuration file for Topup
%__________________________________________________________________________

% Nicole Labra Avila
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

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
% savepwd      = cfg_const;
% savepwd.name = 'Current directory';
% savepwd.tag  = 'savepwd';
% savepwd.val  = {1};
% savepwd.help = {[...
% 'All created files (deformation fields and unwarped images) are ',...
% 'written to the current directory.']};
% 
% saveas       = entry('Save as','ofname','s',[0 Inf]);
% saveas.help  = {[...
% 'Save the result as a VDM file.  "vdm5_" will be prepended to the ',...
% 'filename.']};
% 
% saveusr      = files('Output directory','saveusr','dir',[1 1]);
% saveusr.help = {[...
% 'The combined deformation field and the unwarped images are written ' ...
% 'to the specified directory.']};
% 
% savedir      = cfg_choice;
% savedir.name = 'Output destination';
% savedir.tag  = 'savedir';
% savedir.values = {saveusr savepwd};
% savedir.val  = {saveusr};

%--------------------------------------------------------------------------
% vols Volumes
%--------------------------------------------------------------------------
vols         = cfg_files;
vols.tag     = 'vols';
vols.name    = 'Volumes';
vols.filter  = 'image';
vols.ufilter = '.*';
vols.num     = [2 2];
vols.help    = {[...
'Select two images with opposite phase-encoding polarities. The first one ',...
'is blip up/down and the second one is blip down/up.']};

%--------------------------------------------------------------------------
% fwhm values
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.val     = {[8 4 2 1 0.1]};
fwhm.strtype = 'e';
fwhm.num     = [1 5];
fwhm.help    = {[...
'Spatial scales (expressed as FWHM of Gaussian kernels used to smooth ',...
'input images during topup estimation.']};

%--------------------------------------------------------------------------
% Regularisation 
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Regularisation';
reg.val     = {[0 10 100]}; % Default values 
reg.strtype = 'e';
reg.num     = [1 3];
reg.help    = {[...
'Regularisation settings (see spm_field). The three expected values are:',...
'[1] Penalty on absolute values.',...
'[2] Penalty on the "membrane energy".',...
'[3] Penalty on the "bending energy".']};

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
'a vdm file ("vdm5_*.nii")']};

%--------------------------------------------------------------------------
% disp display
%--------------------------------------------------------------------------
disp         = cfg_menu;
disp.tag     = 'display';
disp.name    = 'Display';
disp.help    = {'Display intermediate results.'};
disp.labels  = {'Yes' 'No'};
disp.values  = {1 0};
disp.val     = {0};

[cfg,varargout{1}] = deal({vols,fwhm,reg,outdir,disp});


%==========================================================================
function out = spm_run_topup(cmd, job)

switch lower(cmd)
    case 'run'
        VDM               = spm_topup(job.vols{1},job.vols{2},job.fwhm, ...
                            job.reg,job.outdir{1},struct('display',job.display));
        out.vdmfile       = {VDM.dat.fname};
        
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Voxel displacement map';
        out(1).src_output = substruct('.','vdmfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
