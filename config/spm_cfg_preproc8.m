function preproc = spm_cfg_preproc8
% Configuration file for 'Combined Segmentation and Spatial Normalisation'
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_preproc8.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% vols Volumes
%--------------------------------------------------------------------------
vols         = cfg_files;
vols.tag     = 'vols';
vols.name    = 'Volumes';
vols.help    = {
    'Select scans from this channel for processing.'
    'If multiple channels are used (eg T1 & T2), then the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..).'
    }';
vols.filter  = 'image';
vols.ufilter = '.*';
vols.num     = [1 Inf];
vols.preview = @(f) spm_check_registration(char(f));

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
biasreg.val    = {0.001};

%--------------------------------------------------------------------------
% biasfwhm Bias FWHM
%--------------------------------------------------------------------------
biasfwhm         = cfg_menu;
biasfwhm.tag     = 'biasfwhm';
biasfwhm.name    = 'Bias FWHM';
biasfwhm.help    = {
    'FWHM of Gaussian smoothness of bias.'
    'If your intensity non-uniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity non-uniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities.'
    }';
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
% write Save Bias Corrected
%--------------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Save Bias Corrected';
write.help    = {
    'Option to save a bias corrected version of your images from this channel, or/and the estimated bias field.'
    'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'
    }';
write.labels = {
                'Save Nothing'
                'Save Bias Corrected'
                'Save Bias Field'
                'Save Field and Corrected'
                }';
write.values = {
                [0 0]
                [0 1]
                [1 0]
                [1 1]
                }';
write.val    = {[0 0]};

%--------------------------------------------------------------------------
% channel Channel
%--------------------------------------------------------------------------
channel         = cfg_branch;
channel.tag     = 'channel';
channel.name    = 'Channel';
channel.val     = {vols biasreg biasfwhm write };
channel.help    = {
    'Specify a channel for processing.'
    'If multiple channels are used (eg PD & T2), then the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..). The different channels can be treated differently in terms of inhomogeneity correction etc. You may wish to correct some channels and save the corrected images, whereas you may wish not to do this for other channels.'
    }';

%--------------------------------------------------------------------------
% data Data
%--------------------------------------------------------------------------
data         = cfg_repeat;
data.tag     = 'data';
data.name    = 'Data';
data.val     = {channel };
data.help    = {
    'Specify the number of different channels (for multi-spectral classification).'
    'If you have scans of different contrasts for each of the subjects, then it is possible to combine the information from them in order to improve the segmentation accuracy. Note that only the first channel of data is used for the initial affine registration with the tissue probability maps.'
    }';
data.values  = {channel };
data.num     = [1 Inf];

%--------------------------------------------------------------------------
% tpm Tissue probability map
%--------------------------------------------------------------------------
tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability map';
tpm.help    = {
    'Select the tissue probability image for this class.'
    'These should be maps of eg grey matter, white matter or cerebro-spinal fluid probability. A nonlinear deformation field is estimated that best overlays the tissue probability maps on the individual subjects'' image.'
    ''
    'Rather than assuming stationary prior probabilities based upon mixing proportions, additional information is used, based on other subjects'' brain images.  Priors are usually generated by registering a large number of subjects together, assigning voxels to different tissue types and averaging tissue classes over subjects. Three tissue classes are used: grey matter, white matter and cerebro-spinal fluid. A fourth class is also used, which is simply one minus the sum of the first three. These maps give the prior probability of any voxel in a registered image being of any of the tissue classes - irrespective of its intensity.'
    ''
    'The model is refined further by allowing the tissue probability maps to be deformed according to a set of estimated parameters. This allows spatial normalisation and segmentation to be combined into the same model.'
    }';
tpm.filter  = 'image';
tpm.ufilter = '.*';
tpm.num     = [1 1];
tpm.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% ngaus Num. Gaussians
%--------------------------------------------------------------------------
ngaus         = cfg_menu;
ngaus.tag     = 'ngaus';
ngaus.name    = 'Num. Gaussians';
ngaus.help    = {
    'The number of Gaussians used to represent the intensity distribution for each tissue class can be greater than one.'
    'In other words, a tissue probability map may be shared by several clusters. The assumption of a single Gaussian distribution for each class does not hold for a number of reasons. In particular, a voxel may not be purely of one tissue type, and instead contain signal from a number of different tissues (partial volume effects). Some partial volume voxels could fall at the interface between different classes, or they may fall in the middle of structures such as the thalamus, which may be considered as being either grey or white matter. Various other image segmentation approaches use additional clusters to model such partial volume effects. These generally assume that a pure tissue class has a Gaussian intensity distribution, whereas intensity distributions for partial volume voxels are broader, falling between the intensities of the pure classes. Unlike these partial volume segmentation approaches, the model adopted here simply assumes that the intensity distribution of each class may not be Gaussian, and assigns belonging probabilities according to these non-Gaussian distributions. Typical numbers of Gaussians could be two for grey matter, two for white matter, two for CSF, three for bone, four for other soft tissues and two for air (background).'
    'Note that if any of the Num. Gaussians is set to non-parametric, then a non-parametric approach will be used to model the tissue intensities. This may work for some images (eg CT), but not others - and it has not been optimised for multi-channel data. Note that it is likely to be especially problematic for images with poorly behaved intensity histograms due to aliasing effects that arise from having discrete values on the images.'
    }';
ngaus.labels = {
                '1'
                '2'
                '3'
                '4'
                '5'
                '6'
                '7'
                '8'
                'Nonparametric'
                }';
ngaus.values = {
                1
                2
                3
                4
                5
                6
                7
                8
                Inf
                }';
ngaus.val    = {Inf};

%--------------------------------------------------------------------------
% native Native Tissue
%--------------------------------------------------------------------------
native         = cfg_menu;
native.tag     = 'native';
native.name    = 'Native Tissue';
native.help    = {'The native space option allows you to produce a tissue class image (c*) that is in alignment with the original/* (see Figure \ref{seg1})*/. It can also be used for ``importing'''' into a form that can be used with the Dartel toolbox (rc*).'};
native.labels = {
    'None'
    'Native Space'
    'Dartel Imported'
    'Native + Dartel Imported'
    }';
native.values = {
                 [0 0]
                 [1 0]
                 [0 1]
                 [1 1]
                 }';
native.val    = {[1 0]};

%--------------------------------------------------------------------------
% warped Warped Tissue
%--------------------------------------------------------------------------
warped         = cfg_menu;
warped.tag     = 'warped';
warped.name    = 'Warped Tissue';
warped.help    = {
    'You can produce spatially normalised versions of the tissue class - both with (mwc*) and without (wc*) modulation (see below). These can be used for voxel-based morphometry. All you need to do is smooth them and do the stats.'
    ''
    '``Modulation'''' is to compensate for the effect of spatial normalisation.  When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images.  For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labelled grey matter.  In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping.  If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.  Actually, in this version of SPM the warped data are not scaled by the Jacobian determinants when generating the "modulated" data.  Instead, the original voxels are projected into their new location in the warped images.  This exactly preserves the tissue count, but has the effect of introducing aliasing artifacts - especially if the original data are at a lower resolution than the warped images.  Smoothing should reduce this artifact though.'
    'Note also that the "unmodulated" data are generated slightly differently in this version of SPM. In this version, the projected data are corrected using a kind of smoothing procedure. This is not done exactly as it should be done (to save computational time), but it does a reasonable job. It also has the effect of extrapolating the warped tissue class images beyond the range of the original data.  This extrapolation is not perfect, as it is only an estimate, but it may still be a good thing to do.'
    }';
warped.labels = {
                 'None'
                 'Modulated'
                 'Unmodulated'
                 'Modulated + Unmodulated'
                 }';
warped.values = {
                 [0 0]
                 [0 1]
                 [1 0]
                 [1 1]
                 }';
warped.val    = {[0 0]};

%--------------------------------------------------------------------------
% tissue Tissue
%--------------------------------------------------------------------------
tissue         = cfg_branch;
tissue.tag     = 'tissue';
tissue.name    = 'Tissue';
tissue.val     = {tpm ngaus native warped};
tissue.help    = {'A number of options are available for each of the tissues.  You may wish to save images of some tissues, but not others. If planning to use Dartel, then make sure you generate ``imported'''' tissue class images of grey and white matter (and possibly others).  Different numbers of Gaussians may be needed to model the intensity distributions of the various tissues.'};

%--------------------------------------------------------------------------
% tissues Tissues
%--------------------------------------------------------------------------
tissues         = cfg_repeat;
tissues.tag     = 'tissues';
tissues.name    = 'Tissues';
tissues.values  = {tissue};
tissues.num     = [0 Inf];

tissues.val     = {tissue tissue tissue tissue tissue tissue};
tpm_nam = fullfile(spm('dir'),'tpm','TPM.nii');
ngaus   = [1 1 2 3 4 2];
% Change to: ngaus   = [2 2 2 3 4 2];
nval    = {[1 0],[1 0],[1 0],[1 0],[1 0],[0 0]};
for k=1:numel(ngaus),
    tissue.val{1}.val = {{[tpm_nam ',' num2str(k)]}};
    tissue.val{2}.val = {ngaus(k)};
    tissue.val{3}.val = {nval{k}};
    tissues.val{k}    = tissue;
end

tissues.help    = {'The data for each subject are classified into a number of different tissue types.  The tissue types are defined according to tissue probability maps, which define the prior probability of finding a tissue type at a particular location. Typically, the order of tissues is grey matter, white matter, CSF, bone, soft tissue and air/background (if using tpm/TPM.nii).'};
tissues.preview = @(f) spm_check_registration(char([f.tpm]));

%--------------------------------------------------------------------------
% mrf MRF Parameter
%--------------------------------------------------------------------------
mrf         = cfg_entry;
mrf.tag     = 'mrf';
mrf.name    = 'MRF Parameter';
mrf.help    = {'When tissue class images are written out, a few iterations of a simple Markov Random Field (MRF) cleanup procedure are run.  This parameter controls the strength of the MRF. Setting the value to zero will disable the cleanup.'};
mrf.strtype = 'r';
mrf.num     = [1 1];
mrf.val     = {1};

%--------------------------------------------------------------------------
% cleanup Clean up any partitions
%--------------------------------------------------------------------------
cleanup         = cfg_menu;
cleanup.tag     = 'cleanup';
cleanup.name    = 'Clean Up';
cleanup.help    = {
    'This uses a crude routine for extracting the brain from segmented images.'
    'It begins by taking the white matter, and eroding it a couple of times to get rid of any odd voxels.  The algorithm continues on to do conditional dilations for several iterations, where the condition is based upon gray or white matter being present.This identified region is then used to clean up the grey and white matter partitions.  Note that the fluid class will also be cleaned, such that aqueous and vitreous humour in the eyeballs, as well as other assorted fluid regions (except CSF) will be removed.'
    ''
    'If you find pieces of brain being chopped out in your data, then you may wish to disable or tone down the cleanup procedure. Note that the procedure uses a number of assumptions about what each tissue class refers to.  If a different set of tissue priors are used, then this routine should be disabled.'
    }';
cleanup.labels = {
                  'Dont do cleanup'
                  'Light Clean'
                  'Thorough Clean'
}';
cleanup.values = {0 1 2};
cleanup.val    = {1};

%--------------------------------------------------------------------------
% reg Warping Regularisation
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Warping Regularisation';
reg.help    = {...
    'Registration involves simultaneously minimising two terms.  One of these is a measure of similarity between the images (mean-squared difference in the current situation), whereas the other is a measure of the roughness of the deformations.  This measure of roughness involves the sum of the following terms:',...
    '* Absolute displacements need to be penalised by a tiny amount.  The first element encodes the amount of penalty on these.  Ideally, absolute displacements should not be penalised, but it is necessary for technical reasons.',...
    '* The `membrane energy'' of the deformation is penalised (2nd element), usually by a relatively small amount. This penalises the sum of squares of the derivatives of the velocity field (ie the sum of squares of the elements of the Jacobian tensors).',...
    '* The `bending energy'' is penalised (3rd element). This penalises the sum of squares of the 2nd derivatives of the velocity.',...
    '* Linear elasticity regularisation is also included (4th and 5th elements).  The first parameter (mu) is similar to that for linear elasticity, except it penalises the sum of squares of the Jacobian tensors after they have been made symmetric (by averaging with the transpose).  This term essentially penalises length changes, without penalising rotations.',...
    '* The final term also relates to linear elasticity, and is the weight that denotes how much to penalise changes to the divergence of the velocities (lambda).  This divergence is a measure of the rate of volumetric expansion or contraction.',...
    'The amount of regularisation determines the tradeoff between the terms. More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations.'};
reg.strtype = 'r';
reg.num     = [1  5];
reg.val     = {[0 0.001 0.5 0.05 0.2]};

%--------------------------------------------------------------------------
% affreg Affine Regularisation
%--------------------------------------------------------------------------
affreg         = cfg_menu;
affreg.tag     = 'affreg';
affreg.name    = 'Affine Regularisation';
affreg.help    = {
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
% fwhm Smoothness
%--------------------------------------------------------------------------
smo         = cfg_entry;
smo.tag     = 'fwhm';
smo.name    = 'Smoothness';
smo.help    = {
    'This is used to derive a fudge factor to account for correlations between neighbouring voxels.'
    'Smoother data have more spatial correlations, rendering the assumptions of the model inaccurate.'
    'For PET or SPECT, set this value to about 5 mm, or more if the images have smoother noise.  For MRI, you can usually use a value of 0 mm.'
    }';
smo.strtype = 'r';
smo.num     = [1  1];
smo.val     = {0};

%--------------------------------------------------------------------------
% samp Sampling distance
%--------------------------------------------------------------------------
samp         = cfg_entry;
samp.tag     = 'samp';
samp.name    = 'Sampling distance';
samp.help    = {
    'This encodes the approximate distance between sampled points when estimating the model parameters.'
    'Smaller values use more of the data, but the procedure is slower and needs more memory. Determining the ``best'''' setting involves a compromise between speed and accuracy.'
    }';
samp.strtype = 'r';
samp.num     = [1  1];
samp.val     = {3};

%--------------------------------------------------------------------------
% write Deformation Fields
%--------------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Deformation Fields';
write.help    = {
    'Deformation fields can be saved to disk, and used by the Deformations Utility.'
    'For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'
    };
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
% warp Warping
%--------------------------------------------------------------------------
warp         = cfg_branch;
warp.tag     = 'warp';
warp.name    = 'Warping & MRF';
warp.val     = {mrf cleanup reg affreg smo samp write};
warp.help    = {
    'A number of warping options.'
    'The main one that you could consider changing is the one for specifying whether deformation fields or inverse deformation fields should be generated.'
    }';

%--------------------------------------------------------------------------
% preproc Segment
%--------------------------------------------------------------------------
preproc      = cfg_exbranch;
preproc.tag  = 'preproc';
preproc.name = 'Segment';
preproc.val  = {data tissues warp};
preproc.help = {
    'Segmentation, bias correction and spatially normalisation - all in the same model.'
    ''
    'This procedure is an extension of the old unified segmentation algorithm (and was known as "New Segment" in SPM8). The algorithm is essentially the same as that described in the Unified Segmentation paper /* \cite{ashburner05}*/, except for (i) a slightly different treatment of the mixing proportions, (ii) the use of an improved registration model, (iii) the ability to use multi-spectral data, (iv) an extended set of tissue probability maps, which allows a different treatment of voxels outside the brain. Some of the options in the toolbox do not yet work, and it has not yet been seamlessly integrated into the SPM8 software.  Also, the extended tissue probability maps need further refinement. The current versions were crudely generated (by JA) using data that was kindly provided by Cynthia Jongen of the Imaging Sciences Institute at Utrecht, NL.'
    ''
    'This function segments, bias corrects and spatially normalises - all in the same model/* \cite{ashburner05}*/. Many investigators use tools within older versions of SPM for a technique that has become known as "optimised" voxel-based morphometry (VBM). VBM performs region-wise volumetric comparisons among populations of subjects. It requires the images to be spatially normalised, segmented into different tissue classes, and smoothed, prior to performing statistical tests/* \cite{wright_vbm,am_vbmreview,ashburner00b,john_should}*/. The "optimised" pre-processing strategy involved spatially normalising subjects'' brain images to a standard space, by matching grey matter in these images, to a grey matter reference.  The historical motivation behind this approach was to reduce the confounding effects of non-brain (e.g. scalp) structural variability on the registration. Tissue classification in older versions of SPM required the images to be registered with tissue probability maps. After registration, these maps represented the prior probability of different tissue classes being found at each location in an image.  Bayes rule can then be used to combine these priors with tissue type probabilities derived from voxel intensities, to provide the posterior probability.'
    ''
    'This procedure was inherently circular, because the registration required an initial tissue classification, and the tissue classification requires an initial registration.  This circularity is resolved here by combining both components into a single generative model. This model also includes parameters that account for image intensity non-uniformity. Estimating the model parameters (for a maximum a posteriori solution) involves alternating among classification, bias correction and registration steps. This approach provides better results than simple serial applications of each component.'
    }';
preproc.prog = @spm_local_preproc_run;
preproc.vout = @vout;


%==========================================================================
function varargout = spm_local_preproc_run(job)
varargout{1} = spm_preproc_run(job);


%==========================================================================
function dep = vout(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

cdep = cfg_dep;
cdep(end).sname      = 'Seg Params';
cdep(end).src_output = substruct('.','param','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

for i=1:numel(job.channel)
    if job.channel(i).write(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('Bias Field (%d)',i);
        cdep(end).src_output = substruct('.','channel','()',{i},'.','biasfield','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.channel(i).write(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('Bias Corrected (%d)',i);
        cdep(end).src_output = substruct('.','channel','()',{i},'.','biascorr','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

for i=1:numel(job.tissue)
    if job.tissue(i).native(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('c%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','c','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).native(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('rc%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','rc','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('wc%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','wc','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('mwc%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','mwc','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

if job.warp.write(1)
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Inverse Deformations';
    cdep(end).src_output = substruct('.','invdef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

if job.warp.write(2)
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Forward Deformations';
    cdep(end).src_output = substruct('.','fordef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

dep = cdep;
