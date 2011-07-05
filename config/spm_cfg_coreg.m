function coreg = spm_cfg_coreg
% SPM Configuration file for Coregister
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_coreg.m 4380 2011-07-05 11:27:12Z volkmar $

% ---------------------------------------------------------------------
% ref Reference Image
% ---------------------------------------------------------------------
ref         = cfg_files;
ref.tag     = 'ref';
ref.name    = 'Reference Image';
ref.help    = {'This is the image that is assumed to remain stationary (sometimes known as the target or template image), while the source image is moved to match it.'};
ref.filter = 'image';
ref.ufilter = '.*';
ref.num     = [1 1];
% ---------------------------------------------------------------------
% source Source Image
% ---------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'This is the image that is jiggled about to best match the reference.'};
source.filter = 'image';
source.ufilter = '.*';
source.num     = [1 1];
% ---------------------------------------------------------------------
% other Other Images
% ---------------------------------------------------------------------
other         = cfg_files;
other.tag     = 'other';
other.name    = 'Other Images';
other.val     = {{''}};
other.help    = {'These are any images that need to remain in alignment with the source image.'};
other.filter = 'image';
other.ufilter = '.*';
other.num     = [0 Inf];
% ---------------------------------------------------------------------
% cost_fun Objective Function
% ---------------------------------------------------------------------
cost_fun         = cfg_menu;
cost_fun.tag     = 'cost_fun';
cost_fun.name    = 'Objective Function';
cost_fun.help    = {'Registration involves finding parameters that either maximise or minimise some objective function. For inter-modal registration, use Mutual Information/* \cite{collignon95,wells96}*/, Normalised Mutual Information/* \cite{studholme99}*/, or Entropy Correlation Coefficient/* \cite{maes97}*/.For within modality, you could also use Normalised Cross Correlation.'};
cost_fun.labels = {
                   'Mutual Information'
                   'Normalised Mutual Information'
                   'Entropy Correlation Coefficient'
                   'Normalised Cross Correlation'
}';
cost_fun.values = {
                   'mi'
                   'nmi'
                   'ecc'
                   'ncc'
}';
cost_fun.def    = @(val)spm_get_defaults('coreg.estimate.cost_fun', val{:});
% ---------------------------------------------------------------------
% sep Separation
% ---------------------------------------------------------------------
sep         = cfg_entry;
sep.tag     = 'sep';
sep.name    = 'Separation';
sep.help    = {'The average distance between sampled points (in mm).  Can be a vector to allow a coarse registration followed by increasingly fine ones.'};
sep.strtype = 'e';
sep.num     = [1 Inf];
sep.def     = @(val)spm_get_defaults('coreg.estimate.sep', val{:});
% ---------------------------------------------------------------------
% tol Tolerances
% ---------------------------------------------------------------------
tol         = cfg_entry;
tol.tag     = 'tol';
tol.name    = 'Tolerances';
tol.help    = {'The accuracy for each parameter.  Iterations stop when differences between successive estimates are less than the required tolerance.'};
tol.strtype = 'e';
tol.num     = [1 12];
tol.def     = @(val)spm_get_defaults('coreg.estimate.tol', val{:});
% ---------------------------------------------------------------------
% fwhm Histogram Smoothing
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Histogram Smoothing';
fwhm.help    = {'Gaussian smoothing to apply to the 256x256 joint histogram. Other information theoretic coregistration methods use fewer bins, but Gaussian smoothing seems to be more elegant.'};
fwhm.strtype = 'e';
fwhm.num     = [1 2];
fwhm.def     = @(val)spm_get_defaults('coreg.estimate.fwhm', val{:});
% ---------------------------------------------------------------------
% eoptions Estimation Options
% ---------------------------------------------------------------------
eoptions         = cfg_branch;
eoptions.tag     = 'eoptions';
eoptions.name    = 'Estimation Options';
eoptions.val     = {cost_fun sep tol fwhm };
eoptions.help    = {'Various registration options, which are passed to the Powell optimisation algorithm/* \cite{press92}*/.'};
% ---------------------------------------------------------------------
% estimate Coreg: Estimate
% ---------------------------------------------------------------------
estimate         = cfg_exbranch;
estimate.tag     = 'estimate';
estimate.name    = 'Coregister: Estimate';
estimate.val     = {ref source other eoptions };
estimate.help    = {
                    'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. The original interpolation method described in this paper has been changed in order to give a smoother cost function.  The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
                    ''
                    'At the end of coregistration, the voxel-to-voxel affine transformation matrix is displayed, along with the histograms for the images in the original orientations, and the final orientations.  The registered images are displayed at the bottom.'
                    ''
                    'Registration parameters are stored in the headers of the "source" and the "other" images.'
}';
estimate.prog = @spm_run_coreg_estimate;
estimate.vout = @vout_estimate;
% ---------------------------------------------------------------------
% ref Image Defining Space
% ---------------------------------------------------------------------
refwrite         = cfg_files;
refwrite.tag     = 'ref';
refwrite.name    = 'Image Defining Space';
refwrite.help    = {'This is analogous to the reference image.  Images are resliced to match this image (providing they have been coregistered first).'};
refwrite.filter = 'image';
refwrite.ufilter = '.*';
refwrite.num     = [1 1];
% ---------------------------------------------------------------------
% source Images to Reslice
% ---------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Images to Reslice';
source.help    = {'These images are resliced to the same dimensions, voxel sizes, orientation etc as the space defining image.'};
source.filter = 'image';
source.ufilter = '.*';
source.num     = [1 Inf];
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Bilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation/* \cite{thevenaz00a,unser93a,unser93b}*/ is slower because it uses more neighbours.'};
interp.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-Spline'
                 '3rd Degree B-Spline'
                 '4th Degree B-Spline'
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values = {0 1 2 3 4 5 6 7};
interp.def     = @(val)spm_get_defaults('coreg.write.interp', val{:});
% ---------------------------------------------------------------------
% wrap Wrapping
% ---------------------------------------------------------------------
wrap         = cfg_menu;
wrap.tag     = 'wrap';
wrap.name    = 'Wrapping';
wrap.help    = {
                'These are typically:'
                '    No wrapping - for PET or images that have already                  been spatially transformed.'
                '    Wrap in  Y  - for (un-resliced) MRI where phase encoding                  is in the Y direction (voxel space).'
}';
wrap.labels = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrap.values = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrap.def     = @(val)spm_get_defaults('coreg.write.wrap', val{:});
% ---------------------------------------------------------------------
% mask Masking
% ---------------------------------------------------------------------
mask         = cfg_menu;
mask.tag     = 'mask';
mask.name    = 'Masking';
mask.help    = {'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};
mask.labels = {
               'Mask images'
               'Dont mask images'
}';
mask.values = {1 0};
mask.def     = @(val)spm_get_defaults('coreg.write.mask', val{:});
% ---------------------------------------------------------------------
% prefix Filename Prefix
% ---------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''r''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('coreg.write.prefix', val{:});
% ---------------------------------------------------------------------
% roptions Reslice Options
% ---------------------------------------------------------------------
roptions         = cfg_branch;
roptions.tag     = 'roptions';
roptions.name    = 'Reslice Options';
roptions.val     = {interp wrap mask prefix };
roptions.help    = {'Various reslicing options.'};
% ---------------------------------------------------------------------
% write Coreg: Reslice
% ---------------------------------------------------------------------
write         = cfg_exbranch;
write.tag     = 'write';
write.name    = 'Coregister: Reslice';
write.val     = {refwrite source roptions };
write.help    = {'Reslice images to match voxel-for-voxel with an image defining some space. The resliced images are named the same as the originals except that they are prefixed by ''r''.'};
write.prog = @spm_run_coreg_reslice;
write.vout = @vout_reslice;
% ---------------------------------------------------------------------
% source Source Image
% ---------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'This is the image that is jiggled about to best match the reference.'};
source.filter = 'image';
source.ufilter = '.*';
source.num     = [1 1];
% ---------------------------------------------------------------------
% estwrite Coreg: Estimate & Reslice
% ---------------------------------------------------------------------
estwrite         = cfg_exbranch;
estwrite.tag     = 'estwrite';
estwrite.name    = 'Coregister: Estimate & Reslice';
estwrite.val     = {ref source other eoptions roptions };
estwrite.help    = {
                    'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. The original interpolation method described in this paper has been changed in order to give a smoother cost function.  The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
                    ''
                    'At the end of coregistration, the voxel-to-voxel affine transformation matrix is displayed, along with the histograms for the images in the original orientations, and the final orientations.  The registered images are displayed at the bottom.'
                    ''
                    'Registration parameters are stored in the headers of the "source" and the "other" images. These images are also resliced to match the source image voxel-for-voxel. The resliced images are named the same as the originals except that they are prefixed by ''r''.'
}';
estwrite.prog = @spm_run_coreg_estwrite;
estwrite.vout = @vout_estwrite;
% ---------------------------------------------------------------------
% coreg Coreg
% ---------------------------------------------------------------------
coreg         = cfg_choice;
coreg.tag     = 'coreg';
coreg.name    = 'Coregister';
coreg.help    = {
                 'Within-subject registration using a rigid-body model. A rigid-body transformation (in 3D) can be parameterised by three translations and three rotations about the different axes.'
                 ''
                 'You get the options of estimating the transformation, reslicing images according to some rigid-body transformations, or estimating and applying rigid-body transformations.'
}';
coreg.values  = {estimate write estwrite };
%coreg.num     = [1 Inf];
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout_estimate(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Coregistered Images';
dep(1).src_output = substruct('.','cfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'Coregistration Matrix';
dep(2).src_output = substruct('.','M');
dep(2).tgt_spec   = cfg_findspec({{'strtype','r'}});
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout_reslice(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Resliced Images';
dep(1).src_output = substruct('.','rfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout_estwrite(job)
depe = vout_estimate(job);
depc = vout_reslice(job);
dep = [depe depc];
