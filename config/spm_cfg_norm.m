function normalise = spm_cfg_norm
% SPM Configuration file for Spatial Normalisation
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_norm.m 6772 2016-04-19 10:21:41Z john $


%--------------------------------------------------------------------------
% biasreg Bias regularisation
%--------------------------------------------------------------------------
biasreg         = cfg_menu;
biasreg.tag     = 'biasreg';
biasreg.name    = 'Bias regularisation';
biasreg.help    = {
                   'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   ''
                   'An important issue relates to the distinction between intensity variations that arise because of bias artifact due to the physics of MR scanning, and those that arise due to different tissue properties.  The objective is to model the latter by different tissue classes, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity non-uniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of non-uniformity.'
                   'Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
                   }';
biasreg.labels = {
                  'no regularisation (0)'
                  'extremely light regularisation (0.00001)'
                  'very light regularisation (0.0001)'
                  'light regularisation (0.001)'
                  'medium regularisation (0.01)'
                  'heavy regularisation (0.1)'
                  'very heavy regularisation (1)'
                  'extremely heavy regularisation (10)'
                  }';
biasreg.values = {
                  0
                  1e-05
                  0.0001
                  0.001
                  0.01
                  0.1
                  1
                  10
                  }';
biasreg.val    = {0.0001};

%--------------------------------------------------------------------------
% biasfwhm Bias FWHM
%--------------------------------------------------------------------------
biasfwhm        = cfg_menu;
biasfwhm.tag    = 'biasfwhm';
biasfwhm.name   = 'Bias FWHM';
biasfwhm.help   = {'FWHM of Gaussian smoothness of bias. If your intensity non-uniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity non-uniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities.'};
biasfwhm.labels = {
                   '30mm cutoff'
                   '40mm cutoff'
                   '50mm cutoff'
                   '60mm cutoff'
                   '70mm cutoff'
                   '80mm cutoff'
                   '90mm cutoff'
                   '100mm cutoff'
                   '110mm cutoff'
                   '120mm cutoff'
                   '130mm cutoff'
                   '140mm cutoff'
                   '150mm cutoff'
                   'No correction'
                   }';
biasfwhm.values = {
                   30
                   40
                   50
                   60
                   70
                   80
                   90
                   100
                   110
                   120
                   130
                   140
                   150
                   Inf
                   }';
biasfwhm.val    = {60};

%--------------------------------------------------------------------------
% write Save Bias Fields
%--------------------------------------------------------------------------
write        = cfg_menu;
write.tag    = 'write';
write.name   = 'Save Bias Fields';
write.help   = {'This is the option concerns whether to save the estimated bias fields. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
write.labels = {
                'Save Nothing'
                'Save Bias Field'
                }';
write.values = {
                [0 0]
                [1 0]
                }';
write.val    = {[0 0]};

%--------------------------------------------------------------------------
% tpm Tissue probability map
%--------------------------------------------------------------------------
tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability map';
tpm.help    = {
               'Select the tissue probability atlas. These should contain probability maps of all the various tissues found in the image data (such that probabilities are greater than or equal to zero, and they sum to one at each voxel. A nonlinear deformation field is estimated that best overlays the atlas on the individual subjects'' image.'
               }';
tpm.filter  = 'nifti';
tpm.ufilter = '.*';
tpm.num     = [1 1];
tpm.val     = {{fullfile(spm('dir'),'tpm','TPM.nii')}};
tpm.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% reg Warping Regularisation
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Warping Regularisation';
reg.help    = {'The objective function for registering the tissue probability maps to the image to process, involves minimising the sum of two terms. One term gives a function of how probable the data is given the warping parameters. The other is a function of how probable the parameters are, and provides a penalty for unlikely deformations. Smoother deformations are deemed to be more probable. The amount of regularisation determines the tradeoff between the terms. Pick a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude). More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations. '};
reg.strtype = 'r';
reg.num     = [1  5];
reg.val     = {[0 0.001 0.5 0.05 0.2]};
%%Eventually these values should be decreased (eg) to:
%reg.val     = {[0 0.0001 0.05 0.005 0.005]};

%--------------------------------------------------------------------------
% affreg Affine Regularisation
%--------------------------------------------------------------------------
affreg        = cfg_menu;
affreg.tag    = 'affreg';
affreg.name   = 'Affine Regularisation';
affreg.help   = {
                  'The procedure is a local optimisation, so it needs reasonable initial starting estimates. Images should be placed in approximate alignment using the Display function of SPM before beginning. A Mutual Information affine registration with the tissue probability maps (D''Agostino et al, 2004) is used to achieve approximate alignment. Note that this step does not include any model for intensity non-uniformity. This means that if the procedure is to be initialised with the affine registration, then the data should not be too corrupted with this artifact.If there is a lot of intensity non-uniformity, then manually position your image in order to achieve closer starting estimates, and turn off the affine registration.'
                  ''
                  'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). For example, if registering to an image in ICBM/MNI space, then choose this option.  If registering to a template that is close in size, then select the appropriate option for this.'
                  }';
affreg.labels = {
                 'No Affine Registration'
                 'ICBM space template - European brains'
                 'ICBM space template - East Asian brains'
                 'Average sized template'
                 'No regularisation'
                 }';
affreg.values = {
                 ''
                 'mni'
                 'eastern'
                 'subj'
                 'none'
                 }';
affreg.val    = {'mni'};

%--------------------------------------------------------------------------
% samp Sampling distance
%--------------------------------------------------------------------------
samp         = cfg_entry;
samp.tag     = 'samp';
samp.name    = 'Sampling distance';
samp.help    = {'This encodes the approximate distance between sampled points when estimating the model parameters. Smaller values use more of the data, but the procedure is slower and needs more memory. Determining the ``best'''' setting involves a compromise between speed and accuracy.'};
samp.strtype = 'r';
samp.num     = [1  1];
samp.val     = {3};

%--------------------------------------------------------------------------
% smo Smoothness
%--------------------------------------------------------------------------
smo         = cfg_entry;
smo.tag     = 'fwhm';
smo.name    = 'Smoothness';
smo.help    = {'For PET or SPECT, set this value to about 5 mm, or more if the images have smoother noise.  For MRI, you can usually use a value of 0 mm.  This is used to derive a fudge factor to account for correlations between neighbouring voxels.  Smoother data have more spatial correlations, rendering the assumptions of the model inaccurate.'};
smo.strtype = 'r';
smo.num     = [1  1];
smo.val     = {0};

%--------------------------------------------------------------------------
% write Deformation Fields
%--------------------------------------------------------------------------
write        = cfg_menu;
write.tag    = 'write';
write.name   = 'Deformation Fields';
write.help   = {'Deformation fields can be saved to disk, and used by the Deformations Utility. For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'};
write.labels = {
                'None'
                'Inverse'
                'Forward'
                'Inverse + Forward'
                }';
write.values = {
                [0 0]
                [1 0]
                [0 1]
                [1 1]
                }';
write.val    = {[0 0]};

%--------------------------------------------------------------------------
% eoptions Estimation Options
%--------------------------------------------------------------------------
eoptions      = cfg_branch;
eoptions.tag  = 'eoptions';
eoptions.name = 'Estimation Options';
eoptions.val  = {biasreg biasfwhm tpm affreg reg smo samp};
eoptions.help = {'Various settings for estimating deformations.'};

%--------------------------------------------------------------------------
% preserve Preserve
%--------------------------------------------------------------------------
preserve        = cfg_menu;
preserve.tag    = 'preserve';
preserve.name   = 'Preserve';
preserve.help   = {
                   'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.'
                   ''
                   'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity.'
}';
preserve.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserve.values = {0 1};
preserve.def    = @(val)spm_get_defaults('normalise.write.preserve', val{:});

%--------------------------------------------------------------------------
% bb Bounding box
%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'};
bb.strtype = 'r';
bb.num     = [2 3];
bb.def     = @(val)spm_get_defaults('normalise.write.bb', val{:});

%--------------------------------------------------------------------------
% vox Voxel sizes
%--------------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel sizes';
vox.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.'};
vox.strtype = 'r';
vox.num     = [1 3];
vox.def     = @(val)spm_get_defaults('normalise.write.vox', val{:});

%--------------------------------------------------------------------------
% interp Interpolation
%--------------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {
                  ['The method by which the images are sampled when ' ...
                  'being written in a different space. ' ...
                  '(Note that Inf or NaN values are treated as zero, ' ...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation' ...
                  '/* \cite{thevenaz00a}*/, especially with higher ' ...
                  'degree splines. Can produce values outside the ' ...
                  'original range (e.g. small negative values from an ' ...
                  'originally all positive image).']
}';
interp.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values = {0 1 2 3 4 5 6 7};
interp.def    = @(val)spm_get_defaults('normalise.write.interp', val{:});

%--------------------------------------------------------------------------
% wrap Wrapping
%--------------------------------------------------------------------------
wrap         = cfg_menu;
wrap.tag     = 'wrap';
wrap.name    = 'Wrapping';
wrap.help    = {
                'These are typically:'
                '    No wrapping: for PET or images that have already been spatially transformed. '
                '    Wrap in  Y: for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'
}';
wrap.labels  = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrap.values  = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrap.def     = @(val)spm_get_defaults('normalise.write.wrap', val{:});

%--------------------------------------------------------------------------
% prefix Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''w''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('normalise.write.prefix', val{:});

%--------------------------------------------------------------------------
% woptions Writing Options
%--------------------------------------------------------------------------
woptions      = cfg_branch;
woptions.tag  = 'woptions';
woptions.name = 'Writing Options';
woptions.val  = { bb vox interp prefix};
woptions.help = {'Various options for writing normalised images.'};

%--------------------------------------------------------------------------
% vol Image to Align
%--------------------------------------------------------------------------
vol         = cfg_files;
vol.tag     = 'vol';
vol.name    = 'Image to Align';
vol.help    = {'The image that the template (atlas) data is warped into alignment with.  The result is a set of warps, which can be applied to this image, or any other image that is in register with it.'};
vol.filter  = 'image';
vol.ufilter = '.*';
vol.num     = [1 1];
vol.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% def Parameter File
%--------------------------------------------------------------------------
def         = cfg_files;
def.tag     = 'def';
def.name    = 'Deformation Field';
def.help    = {[...
    'Deformations can be thought of as vector fields, and represented ',...
    'by three-volume images.  In SPM, deformation fields are saved in ',...
    'NIfTI format, with dimensions xdim x ydim x zdim x 1 x 3. ',...
    'Each voxel contains the x, y and z mm coordinates of where the deformation points.']};
def.filter  = 'nifti';
def.ufilter = 'y_.*\.nii$';
def.num     = [1 1];

%--------------------------------------------------------------------------
% resample Images to Write
%--------------------------------------------------------------------------
resample         = cfg_files;
resample.tag     = 'resample';
resample.name    = 'Images to Write';
resample.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the image used to generate the deformation.'};
resample.filter  = 'image';
resample.ufilter = '.*';
resample.num     = [1 Inf];
resample.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {vol};
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% esubjs Data
%--------------------------------------------------------------------------
esubjs        = cfg_repeat;
esubjs.tag    = 'esubjs';
esubjs.name   = 'Data';
esubjs.help   = {'List of subjects. Images of each subject should be warped differently.'};
esubjs.values = {subj};
esubjs.num    = [1 Inf];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {def resample};
subj.help = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% wsubjs Data
%--------------------------------------------------------------------------
wsubjs        = cfg_repeat;
wsubjs.tag    = 'wsubjs';
wsubjs.name   = 'Data';
wsubjs.help   = {'List of subjects. Images of each subject should be warped differently.'};
wsubjs.values = {subj};
wsubjs.num    = [1 Inf];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {vol resample};
subj.help = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% ewsubjs Data
%--------------------------------------------------------------------------
ewsubjs        = cfg_repeat;
ewsubjs.tag    = 'ewsubjs';
ewsubjs.name   = 'Data';
ewsubjs.help   = {'List of subjects. Images of each subject should be warped differently.'};
ewsubjs.values = {subj};
ewsubjs.num    = [1 Inf];

%--------------------------------------------------------------------------
% est Segment
%--------------------------------------------------------------------------
est      = cfg_exbranch;
est.tag  = 'est';
est.name = 'Normalise: Estimate';
est.val  = {esubjs eoptions};
est.help = {
    'Spatial normalisation is now done via the segmentation routine (which was known as ``New Segment'''' in SPM8).  The algorithm is essentially the same as that described in the Unified Segmentation paper /* \cite{ashburner05}*/, except for (i) a slightly different treatment of the mixing proportions, (ii) the use of an improved registration model, (iii) the ability to use multi-spectral data, (iv) an extended set of tissue probability maps, which allows a different treatment of voxels outside the brain.'
    ''
    'Note that on a 32 bit computer, the most memory that SPM or any other program can use at any time is 4Gbytes (or sometimes only 2Gbytes).  This is because the largest number that can be represented with 32 bits is 4,294,967,295, which limits how much memory may be addressed by any one process.  Out of memory errors may occasionally be experienced when trying to work with large images.  64-bit computers can usually handle such cases.'
        ''
        'If you encounter problems with spatial normalisation, it is advisable to use the Check reg button to see how well aligned the original data are with the MNI-space templates released with SPM.  If mis-alignment is greater than about 3cm and 15 degrees, you could try to manually re-position the images prior to attempting to align them.  This may be done using the Display button.'
           }';
est.prog = @spm_run_norm;
est.vout = @vout_est;

%--------------------------------------------------------------------------
% write Normalise: Write
%--------------------------------------------------------------------------
write      = cfg_exbranch;
write.tag  = 'write';
write.name = 'Normalise: Write';
write.val  = {wsubjs woptions};
write.help = {'Allows previously estimated warps (stored in ``y_''''imagename``_sn.mat'''' files) to be applied to series of images.'};
write.prog = @spm_run_norm;
write.vout = @vout_write;

%--------------------------------------------------------------------------
% estwrite Normalise: Estimate & Write
%--------------------------------------------------------------------------
estwrite      = cfg_exbranch;
estwrite.tag  = 'estwrite';
estwrite.name = 'Normalise: Estimate & Write';
estwrite.val  = {ewsubjs eoptions woptions};
estwrite.help = {'Computes the warp that best aligns the template (atlas) to the individual''s image, inverting it and writing the result to the file `y_''imagename''.nii''. This option also allows the contents of the `y_''imagename''.nii'' files to be applied to a series of images.'
''
'Note that if you encounter problems with spatial normalisation, it is often advisable to use the Check reg button to see how well aligned the original data are with the MNI-space templates released with SPM.  If mis-alignment is greater than about 3cm and 15 degrees, you could try to manually re-position the images.  This may be done using the Display button.'};
estwrite.prog = @spm_run_norm;
estwrite.vout = @vout_estwrite;

%--------------------------------------------------------------------------
% normalise Normalise
%--------------------------------------------------------------------------
normalise        = cfg_choice;
normalise.tag    = 'normalise';
normalise.name   = 'Normalise';
normalise.help   = {...
    'There are two components to spatial normalisation: There is the estimation part, ',...
    'whereby a deformation is estimated by deforming template data to match an ',...
    'individual scan; And there is the actual writing of the spatially normalised ',...
    'images, using the previously estimated deformation.',...
    'This is a vanilla approach to spatial normalisation.  ',...
    'It is not generally recommended for morphometric studies, or other studies of ',...
    'differences among populations. ',...
    'The reason is that the atlas data will differ systematically from the data under study, ',...
    'which is likely to lead to an inherently biased set of findings.' };
normalise.values = {est write estwrite};


%==========================================================================
function dep = vout_est(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Deformation (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','def');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end


%==========================================================================
function dep = vout_write(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Normalised Images (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end


%==========================================================================
function dep = vout_estwrite(job)
depe = vout_est(job);
depw = vout_write(job);
dep = [depe depw];
