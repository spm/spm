function slice2vol = spm_slice2vol_config
% Configuration file for toolbox 'Spatial Tools'
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Spatial')); end

% -------------------------------------------------------------------------
slice2vol      = cfg_exbranch;
slice2vol.tag  = 'slice2vol';
slice2vol.name = 'Slice-to-volume alignment';
slice2vol.val  = @slice2vol_cfg;
slice2vol.help = {...
'Slice-to-volume alignment of fMRI time series.',...
['The rigid-body alignment assumption is a bit limited for fMRI data, ',...
 'which is normally acquired a slice at a time. ',...
 'Motion during the acquisition of a volume is not properly accounted ',...
 'for with volume-based alignment, so this slice-to-volume alignment ',...
 'is an approach to try to deal with this.'],...
['Note that the method is a little slower than conventional motion ',...
 'correction, because it uses more of the data in order to estimate rather ',...
 'more parameters.'],...
['Also note that certain types of motion can cause the slices to move ',...
 'such that some of the data in the resliced volume is missing. ',...
 'This is why smoothing options are included here.'] };
slice2vol.prog = @spm_slice2vol;
slice2vol.vout = @vout_slice2vol;


function varargout = slice2vol_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select fMRI data to correct.'};
images.filter = 'nifti';
images.ufilter = '^.*';
images.num     = [1 Inf];
images.preview = @(f) spm_check_registration(char(f));
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
sd         = cfg_entry;
sd.tag     = 'sd';
sd.name    = 'Motion sd';
sd.help    = {...
['Motion parameters are regularised according to how much motion (mm) ',...
 'might typically be expected from one slice to the next. This is to prevent ',...
 'overfitting because there are many parameters to be estimated.']...
};
sd.strtype = 'r';
sd.num     = [1 1];
sd.val     = {0.001};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
slice_code        = cfg_menu;
slice_code.tag    = 'slice_code';
slice_code.name   = 'Slice Ordering';
slice_code.help   = {[...
'Motion parameter estimates are regularised so that slices adjacent in ',...
'time have to move in a similar way to each other. Slice ordering is therefore ',...
'used to determine temporal adjacency.']};
slice_code.labels = {
                 'Sequential Increasing'
                 'Sequential Decreasing'
                 'Alternating Increasing'
                 'Alternating Decreasing'
                 'Unknown'
}';
slice_code.values = {1 2 3 4 0};
slice_code.val    = {1};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {...
['Full width at half maximum (FWHM) of the isotropic Gaussian smoothing ',...
 'kernel in mm. This is used to smooth the motion corrected images.']...
};
fwhm.strtype = 'r';
fwhm.num     = [1 1];
fwhm.val     = {0};
%--------------------------------------------------------------------------

[cfg,varargout{1}] = deal({images,slice_code,sd,fwhm});


%==========================================================================
function dep = vout_slice2vol(job)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'Mean Image';
dep(1).src_output = substruct('.','rmean');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Resliced Images';
dep(2).src_output = substruct('.','rfiles');
dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
