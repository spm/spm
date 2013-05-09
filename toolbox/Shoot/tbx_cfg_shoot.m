function shoot = tbx_cfg_shoot
% MATLABBATCH Configuration file for toolbox 'Shoot Tools'

% $Id: tbx_cfg_shoot.m 5485 2013-05-09 15:51:24Z john $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Shoot')); end

% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select a set of imported images of the same type to be registered by minimising a measure of difference from the template.'};
images1.filter = 'image';
images1.ufilter = '^r.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the images to be warped together. Multiple sets of images can be simultaneously registered. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.'};
images.values  = {images1};
images.val     = {images1};
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% warp Run Shoot (create Templates)
% ---------------------------------------------------------------------
warp         = cfg_exbranch;
warp.tag     = 'warp';
warp.name    = 'Run Shooting (create Templates)';
warp.val     = {images };
warp.help    = {'Run the geodesic shooting nonlinear image registration procedure. This involves iteratively matching all the selected images to a template generated from their own mean. A series of Template*.nii files are generated, which become increasingly crisp as the registration proceeds.'};
warp.prog = @spm_shoot_template;
warp.vout = @vout_shoot_template;
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select a set of imported images of the same type to be registered by minimising a measure of difference from the template.'};
images1.filter = 'image';
images1.ufilter = '^r.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the images to be warped together. Multiple sets of images can be simultaneously registered. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'templates';
template.name    = 'Templates';
template.help    = {'Select templates. Smoother templates should be used for the early iterations. Note that the template should be a 4D file, with the 4th dimension equal to the number of sets of images.'};
template.filter = 'nifti';
template.ufilter = '.*';
template.num     = [1 Inf];
% ---------------------------------------------------------------------
% warp1 Run Shooting (existing Templates)
% ---------------------------------------------------------------------
warp1         = cfg_exbranch;
warp1.tag     = 'warp1';
warp1.name    = 'Run Shoot (existing Templates)';
warp1.val     = {images template };
warp1.check   = @check_shoot_template;
warp1.help    = {'Run the Shoot nonlinear image registration procedure to match individual images to pre-existing template data. Start out with smooth templates, and select crisp templates for the later iterations.'};
warp1.prog = @spm_shoot_warp;
warp1.vout = @vout_shoot_warp;
% ---------------------------------------------------------------------
% velocities Velocity fields
% ---------------------------------------------------------------------
velocities         = cfg_files;
velocities.tag     = 'velocities';
velocities.name    = 'Velocity fields';
velocities.help    = {'The velocity fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
velocities.filter = 'nifti';
velocities.ufilter = '^v_.*';
velocities.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are velocity fields, such that each velocity field warps one image.'};
images1.filter = 'nifti';
images1.ufilter = '.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'The velocity field deformations can be applied to multiple images. At this point, you are choosing how many images each velocity field should be applied to.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% See later
% ---------------------------------------------------------------------
many_subj = cfg_branch;
many_subj.tag = 'subjs';
many_subj.name = 'Many Subjects';
many_subj.val  = {velocities,images};
many_subj.help = {[...
'Select this option if you have many subjects to spatially normalise, ',...
'but there are a small and fixed number of scans for each subject.']};
% ---------------------------------------------------------------------
% jactransf Modulation
% ---------------------------------------------------------------------
jactransf         = cfg_menu;
jactransf.tag     = 'jactransf';
jactransf.name    = 'Modulation';
jactransf.val     = {0};
jactransf.help    = {'This allows the spatially normalised images to be rescaled by the Jacobian determinants of the deformations. Note that the rescaling is only approximate for deformations generated using smaller numbers of time steps.'};
jactransf.labels  = {
                    'Pres. Concentration (No "modulation")'
                    'Pres. Amount ("Modulation")'
}';
jactransf.values  = {0 1};
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.val     = {1};
interp.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:          - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:    - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Do not use B-splines when there is any region of NaN or Inf in the images. '
}';
interp.labels  = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values  = {0 1 2 3 4 5 6 7};
% ---------------------------------------------------------------------
% crt_warped Create Warped
% ---------------------------------------------------------------------
crt_warped         = cfg_exbranch;
crt_warped.tag     = 'crt_warped';
crt_warped.name    = 'Create Warped';
crt_warped.val     = {velocities images jactransf interp };
crt_warped.check   = @check_norm;
crt_warped.help    = {'This allows spatially normalised images to be generated. Note that voxel sizes and bounding boxes can not be adjusted, and that there may be strange effects due to the boundary conditions used by the warping. Also note that the warped images are not in Talairach or MNI space. The coordinate system is that of the average shape and size of the subjects to which Shoot was applied. In order to have MNI-space normalised images, then the Deformations Utility can be used to compose the individual Shoot warps, with a deformation field that matches (e.g.) the Template grey matter generated by Shoot, with one of the grey matter volumes released with SPM.'};
%crt_warped.prog = @spm_shoot_norm;
%crt_warped.vout = @vout_norm;
% ---------------------------------------------------------------------
% velocities Velocity fields
% ---------------------------------------------------------------------
velocities         = cfg_files;
velocities.tag     = 'velocities';
velocities.name    = 'Velocity fields';
velocities.help    = {'The velocity fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
velocities.filter = 'nifti';
velocities.ufilter = '^v_.*';
velocities.num     = [1 Inf];
% ---------------------------------------------------------------------
% jacdet Jacobian determinants
% ---------------------------------------------------------------------
jacdet         = cfg_exbranch;
jacdet.tag     = 'jacdet';
jacdet.name    = 'Jacobian determinants';
jacdet.val     = {velocities};
jacdet.help    = {'Create Jacobian determinant fields from velocities.'};
%jacdet.prog = @spm_shoot_jacobian;
% ---------------------------------------------------------------------
% velocities Velocity fields
% ---------------------------------------------------------------------
velocities         = cfg_files;
velocities.tag     = 'velocities';
velocities.name    = 'Velocity fields';
velocities.help    = {'The velocity fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
velocities.filter = 'nifti';
velocities.ufilter = '^v_.*';
velocities.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the image(s) to be inverse normalised.  These should be in alignment with the template image generated by the warping procedure.'};
images.filter = 'nifti';
images.ufilter = '.*';
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.val{1} = double(1);
interp.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:          - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:    - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Do not use B-splines when there is any region of NaN or Inf in the images.'
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
interp.values  = {0 1 2 3 4 5 6 7};
% ---------------------------------------------------------------------
% crt_iwarped Create Inverse Warped
% ---------------------------------------------------------------------
crt_iwarped         = cfg_exbranch;
crt_iwarped.tag     = 'crt_iwarped';
crt_iwarped.name    = 'Create Inverse Warped';
crt_iwarped.val     = {velocities images interp };
crt_iwarped.help    = {'Create inverse normalised versions of some image(s). The image that is inverse-normalised should be in alignment with the template (generated during the warping procedure). Note that the results have the same dimensions as the ``velocity fields'''', but are mapped to the original images via the affine transformations in their headers.'};
%crt_iwarped.prog = @spm_shoot_invnorm;
%crt_iwarped.vout = @vout_invnorm;
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
velocityfield         = cfg_files;
velocityfield.tag     =  'velocityfield';
velocityfield.name    = 'Velocity Field';
velocityfield.filter  = 'nifti';
velocityfield.ufilter = '^v_.*\.nii$';
velocityfield.num     = [1 1];
velocityfield.help    = {'Shoot velocity field for this subject.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
images        = cfg_files;
images.tag    = 'images';
images.name   = 'Images';
images.filter = 'nifti';
images.num    = [1 Inf];
images.help   = {'Images for this subject to spatially normalise.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {velocityfield,images};
subj.help = {'Subject to be spatially normalized.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
few_subj        = cfg_repeat;
few_subj.tag    = 'few_subj';
few_subj.name   = 'Few Subjects';
few_subj.values = {subj};
few_subj.help   = {[...
'Select this option if there are only a few subjects, each with many or ',...
'a variable number of scans each. You will then need to specify a series of subjects, and the velocity field and images of each of them.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
way        = cfg_choice;
way.tag    = 'data';
way.name   = 'Select according to';
way.values = {few_subj,many_subj};
way.help   = {...
['You may wish to spatially normalise only a few subjects, '...
 'but have many scans per subject (eg for fMRI), '...
 'or you may have lots of subjects, but with a small and fixed number '...
 'of scans for each of them (eg for VBM).  The idea is to chose the way of '...
 'selecting files that is easier.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
template        = cfg_files;
template.tag    = 'template';
template.name   = 'Shoot Template';
template.filter = 'nifti';
template.num    = [0 1];
template.help   = {...
['Select the final Template file generated by Shoot. This will be affine '...
 'registered with a TPM file, such that the resulting spatially normalised '...
 'images are closer aligned to MNI space. Leave empty if you do not wish to '...
 'incorporate a transform to MNI space '...
 '(ie just click ``done'' on the file selector, without selecting any images).']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Gaussian FWHM';
fwhm.val     = {[8 8 8]};
fwhm.strtype = 'e';
fwhm.num     = [1 3];
fwhm.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian blurring kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions. Note that you can also specify [0 0 0], but any ``modulated'' data will show aliasing (see eg Wikipedia), which occurs because of the way the warped images are generated.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
preserve         = cfg_menu;
preserve.tag     = 'preserve';
preserve.name    = 'Preserve';
preserve.help    = {
'Preserve Concentrations: Smoothed spatially normalised images (sw*) represent weighted averages of the signal under the smoothing kernel, approximately preserving the intensities of the original images. This option is currently suggested for eg fMRI.'
''
'Preserve Total: Smoothed and spatially normalised images preserve the total amount of signal from each region in the images (smw*). Areas that are expanded during warping are correspondingly reduced in intensity. This option is suggested for VBM.'
}';
preserve.labels = {
                   'Preserve Concentrations (no "modulation")'
                   'Preserve Amount ("modulation")'
}';
preserve.values = {0 1};
preserve.val    = {0};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
vox          = cfg_entry;
vox.tag      = 'vox';
vox.name     = 'Voxel sizes';
vox.num      = [1 3];
vox.strtype  = 'e';
vox.def      = @(val)spm_get_defaults('defs.vox',val{:});
vox.help     = {[...
'Specify the voxel sizes of the deformation field to be produced. ',...
'Non-finite values will default to the voxel sizes of the template image',...
'that was originally used to estimate the deformation.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
bb           = cfg_entry;
bb.tag       = 'bb';
bb.name      = 'Bounding box';
bb.strtype   = 'e';
bb.num       = [2 3];
bb.def       = @(val)spm_get_defaults('defs.bb',val{:});
bb.help      = {[...
'Specify the bounding box of the deformation field to be produced. ',...
'Non-finite values will default to the bounding box of the template image',...
'that was originally used to estimate the deformation.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
nrm       = cfg_exbranch;
nrm.tag   = 'mni_norm';
nrm.name  = 'Normalise to MNI Space';
nrm.val   = {template,way,vox,bb,preserve,fwhm};
%nrm.prog  = @spm_shoot_norm_fun;
%nrm.vout  = @vout_norm_fun;
%nrm.check = @check_norm_fun;
nrm.help  = {[...
'Normally, Shoot generates warped images that align with the average-shaped template. ',...
'This routine includes an initial affine regisration of the template (the final one ',...
'generated by Shoot), with the TPM data released with SPM.'],[...
'``Smoothed'''' (blurred) spatially normalised images are generated in such a ',...
'way that the original signal is preserved. Normalised images are ',...
'generated by a ``pushing'''' rather than a ``pulling'''' (the usual) procedure. ',...
'Note that trilinear interpolation is used, and no masking is done.  It ',...
'is therefore essential that the images are realigned and resliced ',...
'before they are spatially normalised.  Alternatively, contrast images ',...
'generated from unsmoothed native-space fMRI/PET data can be spatially ',...
'normalised for a 2nd level analysis.'],[...
'Two ``preserve'''' options are provided.  One of them should do the ',...
'equavalent of generating smoothed ``modulated'''' spatially normalised ',...
'images.  The other does the equivalent of smoothing the modulated ',...
'normalised fMRI/PET, and dividing by the smoothed Jacobian determinants.']};
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select tissue class images (one per subject).'};
images1.filter = 'nifti';
images1.ufilter = '^r.*c';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Multiple sets of images are used here. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.  The number of sets of images must be the same as was used to generate the template.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% deformations Deformation fields
% ---------------------------------------------------------------------
deformations         = cfg_files;
deformations.tag     = 'deformations';
deformations.name    = 'Deformation fields';
deformations.help    = {'Select the deformation fields for each subject.'};
deformations.filter = 'nifti';
deformations.ufilter = '^y_.*';
deformations.num     = [1 Inf];
% ---------------------------------------------------------------------
% jacobians Jacobian determinant fields
% ---------------------------------------------------------------------
jacobians         = cfg_files;
jacobians.tag     = 'jacobians';
jacobians.name    = 'Jacobian determinant fields';
jacobians.help    = {'Select the Jacobian determinant fields for each subject.  Residual differences are computed between the warped images and template. These are then scaled by the Jacobian determinants at each point, and spatially smoothed.'};
jacobians.filter = 'nifti';
jacobians.ufilter = '^j_.*';
jacobians.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template';
template.help    = {'Residual differences are computed between the warped images and template. These are then scaled by the Jacobian determinants at each point, and spatially smoothed.'};
template.filter = 'nifti';
template.ufilter = '^Template.*';
template.num     = [0 1];
% ---------------------------------------------------------------------
% fwhm Smoothing
% ---------------------------------------------------------------------
fwhm         = cfg_menu;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Smoothing';
fwhm.val     = {4};
fwhm.help    = {'The scalar momenta can be smoothed with a Gaussian to reduce dimensionality. More smoothing is recommended if there are fewer training images or if more channels of data were used for driving the registration. From preliminary experimants, a value of about 10mm seems to work reasonably well.'};
fwhm.labels  = {
               'None'
               ' 2mm'
               ' 4mm'
               ' 6mm'
               ' 8mm'
               '10mm'
               '12mm'
               '14mm'
               '16mm'
}';
fwhm.values  = {0 2 4 6 8 10 12 14 16};
fwhm.val     = {10};
% ---------------------------------------------------------------------
% scalmom Generate Scalar Momenta
% ---------------------------------------------------------------------
scalmom         = cfg_exbranch;
scalmom.tag     = 'scalmom';
scalmom.name    = 'Generate Scalar Momenta';
scalmom.val     = {template images deformations jacobians fwhm};
scalmom.check   = @check_scalmom;
scalmom.help    = {'Generate spatially smoothed ``scalar momenta'''' /* cite{singh2010multivariate,singh2012genetic} */ in a form suitable for using with pattern recognition. In principle, a Gaussian Process model can be used to determine the optimal (positive) linear combination of kernel matrices.  The idea would be to combine a kernel matrix derived from these, with a kernel derived from the velocity-fields. Such a combined kernel should then encode more relevant information than the individual kernels alone.  The scalar momentum fields that are generated contain a number of volumes equal to the number of sets of ``rc*'''' images used (equal to the number of volumes in the template - 1).  /* See Figures 10 and 11 of \cite{ashburner2011multivariate} for examples of scalar momenta (Jacobian scaled residuals) for simulated data. */'};
scalmom.prog = @spm_shoot_scalmom;
scalmom.vout = @vout_scalmom;
% ---------------------------------------------------------------------
% images Data
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Data';
images.help    = {'Select images to generate dot-products from.'};
images.filter = 'nifti';
images.ufilter = '.*';
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% weight Weighting image
% ---------------------------------------------------------------------
weight         = cfg_files;
weight.tag     = 'weight';
weight.name    = 'Weighting image';
weight.val = {{}};
weight.help    = {'The kernel can be generated so that some voxels contribute to the similarity measures more than others.  This is achieved by supplying a weighting image, which each of the component images are multiplied before the dot-products are computed. This image needs to have the same dimensions as the component images, but orientation information (encoded by matrices in the headers) is ignored. If left empty, then all voxels are weighted equally.'};
weight.filter = 'image';
weight.ufilter = '.*';
weight.num     = [0 1];
% ---------------------------------------------------------------------
% dotprod Dot-product Filename
% ---------------------------------------------------------------------
dotprod         = cfg_entry;
dotprod.tag     = 'dotprod';
dotprod.name    = 'Dot-product Filename';
dotprod.help    = {'Enter a filename for results (it will be prefixed by ``dp_'''' and saved in the current directory).'};
dotprod.strtype = 's';
dotprod.num     = [1 Inf];
% ---------------------------------------------------------------------
% reskern Kernel from Resids
% ---------------------------------------------------------------------
reskern         = cfg_exbranch;
reskern.tag     = 'reskern';
reskern.name    = 'Kernel from Images';
reskern.val     = {images weight dotprod };
reskern.help    = {'Generate a kernel matrix from images. In principle, this same function could be used for generating kernels from any image data (e.g. ``modulated'''' grey matter). If there is prior knowledge about some region providing more predictive information (e.g. the hippocampi for AD), then it is possible to weight the generation of the kernel accordingly. The matrix of dot-products is saved in a variable ``K'''', which can be loaded from the dp_*.mat file. The ``kernel trick'''' can be used to convert these dot-products into distance measures for e.g. radial basis-function approaches.'};
reskern.prog = @spm_dotprods2;
reskern.vout = @vout_dotprod;
% ---------------------------------------------------------------------
% velocities Velocity fields
% ---------------------------------------------------------------------
velocities         = cfg_files;
velocities.tag     = 'velocities';
velocities.name    = 'Velocity fields';
velocities.help    = {'Select the velocity fields for each subject.'};
velocities.filter = 'nifti';
velocities.ufilter = '^v_.*';
velocities.num     = [1 Inf];
% ---------------------------------------------------------------------
% dotprod Dot-product Filename
% ---------------------------------------------------------------------
dotprod         = cfg_entry;
dotprod.tag     = 'dotprod';
dotprod.name    = 'Dot-product Filename';
dotprod.help    = {'Enter a filename for results (it will be prefixed by ``dp_'''' and saved in the current directory.'};
dotprod.strtype = 's';
dotprod.num     = [1 Inf];
% ---------------------------------------------------------------------
% velkern Kernel from Velocities
% ---------------------------------------------------------------------
velkern         = cfg_exbranch;
velkern.tag     = 'velkern';
velkern.name    = 'Kernel from velocities';
velkern.val     = {velocities dotprod};
velkern.help    = {'Generate a kernel from velocity fields. The dot-products are saved in a variable ``K'''' in the resulting dp_*.mat file.'};
velkern.prog = @spm_shoot_kernel;
velkern.vout = @vout_kernel;
% ---------------------------------------------------------------------
% kernfun Kernel Utilities
% ---------------------------------------------------------------------
kernfun         = cfg_choice;
kernfun.tag     = 'kernfun';
kernfun.name    = 'Kernel Utilities';
kernfun.help    = {
                   'Shoot can be used for generating matrices of dot-products for various kernel pattern-recognition procedures.'
                   'The idea of applying pattern-recognition procedures is to obtain a multi-variate characterisation of the anatomical differences among groups of subjects. These characterisations can then be used to separate (eg) healthy individuals from particular patient populations. There is still a great deal of methodological work to be done, so the types of kernel that can be generated here are unlikely to be the definitive ways of proceeding.  They are only just a few ideas that may be worth trying out. The idea is simply to attempt a vaguely principled way to combine generative models with discriminative models (see the ``Pattern Recognition and Machine Learning'''' book by Chris Bishop for more ideas). Better ways (higher predictive accuracy) will eventually emerge.'
                   'Various pattern recognition algorithms are available freely over the Internet. Possible approaches include Support-Vector Machines and Gaussian Process Models. Gaussian Process Models probably give the most accurate probabilistic predictions, and allow kernels generated from different pieces of data to be most easily combined.'
}';
kernfun.values  = {velkern scalmom reskern};
% ---------------------------------------------------------------------
% shoot Shoot Tools
% ---------------------------------------------------------------------
shoot         = cfg_choice;
shoot.tag     = 'shoot';
shoot.name    = 'Shoot Tools';
shoot.help    = {
                  'This toolbox is based around the ``Diffeomorphic Registration using Geodesic Shooting and Gauss-Newton Optimisation'''' paper, which has been submitted to NeuroImage. The idea is to register images by estimating an initial velocity field, which can then be integrated to generate both forward and backward deformations.  Currently, the software only works with images that have isotropic voxels, identical dimensions and which are in approximate alignment with each other. One of the reasons for this is that the approach assumes circulant boundary conditions, which makes modelling global rotations impossible. Because of these limitations, the registration should be based on images that have first been ``imported'''' via the New Segment toolbox.'
                  'The next step is the registration itself, which involves the simultaneous registration of e.g. GM with GM, WM with WM and 1-(GM+WM) with 1-(GM+WM) (when needed, the 1-(GM+WM) class is generated implicitly, so there is no need to include this class yourself). This procedure begins by creating a mean of all the images, which is used as an initial template. Deformations from this template to each of the individual images are computed, and the template is then re-generated by applying the inverses of the deformations to the images and averaging. This procedure is repeated a number of times.'
                  ''
                  'This toolbox should be considered as only a beta (trial) version, and will include a number of (as yet unspecified) extensions in future updates.  Please report any bugs or problems to the SPM mailing list.'
}';
shoot.values  = {warp warp1 kernfun};

%_______________________________________________________________________
%
%_______________________________________________________________________
function chk = check_shoot_template(job)
n1 = numel(job.images);
n2 = numel(job.images{1});
chk = '';
for i=1:n1,
    if numel(job.images{i}) ~= n2,
        chk = 'Incompatible number of images';
        break;
    end;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_shoot_template(job)

d = spm_shoot_defaults;
if isfield(d, 'tname') & ~isempty(deblank(d.tname)),
    for it=0:ceil((numel(d.sched)-1)/6),
        tdep(it+1)            = cfg_dep;
        tdep(it+1).sname      = sprintf('Template (%d)', it);
        tdep(it+1).src_output = substruct('.','template','()',{it+1});
        tdep(it+1).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
else
    tdep = cfg_dep;
    tdep = tdep(false);
end
vdep            = cfg_dep;
vdep.sname      = 'Velocity Fields';
vdep.src_output = substruct('.','vel','()',{':'});
vdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

ydep            = cfg_dep;
ydep.sname      = 'Deformation Fields';
ydep.src_output = substruct('.','def','()',{':'});
ydep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

jdep            = cfg_dep;
jdep.sname      = 'Jacobian Fields';
jdep.src_output = substruct('.','jac','()',{':'});
jdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

dep = [tdep vdep ydep jdep];
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_shoot_warp(job)
dep            = cfg_dep;
dep.sname      = 'Velocity Fields';
dep.src_output = substruct('.','files','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_norm(job)
chk = '';
PU = job.velocities;
PI = job.images;
n1 = numel(PU);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_norm_fun(job)
chk = '';
if isfield(job.data,'subjs')
    PU = job.data.subjs.velocities;
    PI = job.data.subjs.images;
    n1 = numel(PU);
    for i=1:numel(PI),
        if numel(PI{i}) ~= n1,
            chk = 'Incompatible number of images';
            break;
        end
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_norm(job)
if job.jactransf,
    sname = 'Warped Images - Jacobian Transformed';
else
    sname = 'Warped Images';
end
PU    = job.velocities;
PI    = job.images;
for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('%s (Image %d)',sname,m);
    dep(m).src_output = substruct('.','files','()',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
n = numel(PI);
for m=1:numel(PU),
    dep(m+n)            = cfg_dep;
    dep(m+n).sname      = sprintf('%s (Deformation %d)',sname,m);
    dep(m+n).src_output = substruct('.','files','()',{m,':'});
    dep(m+n).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_invnorm(job)
PU    = job.velocities;
PI    = job.images;

for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('Inverse Warped Images (Image %d)',m);
    dep(m).src_output = substruct('.','files','()',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
n = numel(PI);
for m=1:numel(PU),
    dep(m+n)            = cfg_dep;
    dep(m+n).sname      = sprintf('Inverse Warped Images (Deformation %d)',m);
    dep(m+n).src_output = substruct('.','files','()',{m,':'});
    dep(m+n).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_norm_fun(job)
if job.preserve,
    sname = 'MNI Smo. Warped - Amount';
else
    sname = 'MNI Smo. Warped - Concentrations';
end

dep = cfg_dep;

if isfield(job.data,'subj')
    for m=1:numel(job.data.subj)
        dep(m)            = cfg_dep;
        dep(m).sname      = sprintf('%s (Deformation %d)',sname,m);
        dep(m).src_output = substruct('{}',{m},'()',{':'});
        dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

if isfield(job.data,'subjs')
    for m=1:numel(job.data.subjs.images),
        dep(m)            = cfg_dep;
        dep(m).sname      = sprintf('%s (Image %d)',sname,m);
        dep(m).src_output = substruct('()',{':',m});
        dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_scalmom(job)
chk = '';
PY = job.deformations;
PI = job.images;
n1 = numel(PY);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_scalmom(job)
dep            = cfg_dep;
dep.sname      = 'Scalar Momentum Fields';
dep.src_output = substruct('.','scalmom','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_kernel(job)
dep            = cfg_dep;
dep.sname      = 'Velocity Kernel';
dep.src_output = substruct('.','fname','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_dotprod(job)
dep            = cfg_dep;
dep.sname      = 'Image Kernel';
dep.src_output = substruct('.','fname','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
%_______________________________________________________________________

