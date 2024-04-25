function cfg = tbx_cfg_mb
% MATLABBATCH Configuration file for toolbox 'Multi-Brain'
%_____________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

%--------------------------------------------------------------------------
images        = cfg_files;
images.tag    = 'images';
images.name   = 'Scans';
images.filter = 'nifti';
images.num    = [1 Inf];
images.help   = {['Select one NIfTI format scan for each subject. Subjects must be in the same order if there are multiple channels. ' ...
                  'Image dimensions can differ over subjects, but (if there are multiple channels) the scans of each ' ...
                  'subject must all have the same dimensions and orientations.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
inu_reg      = cfg_entry;
inu_reg.tag  = 'inu_reg';
inu_reg.name = 'Regularisation';
inu_reg.strtype = 'e';
inu_reg.num  = [1 1];
inu_reg.val  = {1e4};
inu_reg.help = {['Specify the bending energy penalty on the estimated intensity nonuniformity (INU) ' ...
                'fields (bias fields). Larger values give smoother INU fields.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
inu_co        = cfg_menu;
inu_co.tag    = 'inu_co';
inu_co.name   = 'Cut off';
inu_co.labels = {' 20 mm INU', ' 40 mm INU', ' 60 mm INU', ' 80 mm INU', '100 mm INU', 'Rescale only', 'No correction'};
inu_co.values = {20, 40, 60, 80, 100, Inf, NaN};
inu_co.val    = {40};
inu_co.help   = {['Specify the cutoff (mm) of the intensity nonuniformity (INU) correction (bias correction). ' ...
                 'Larger values use fewer parameters to encode the INU field. ' ...
                 'Note that a global intensity rescaling correction, without INU correction, can also be specified. ' ...
                 'For quantitative images, it may be better not to use any correction.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
inu           = cfg_branch;
inu.tag       = 'inu';
inu.name      = 'Intensity nonuniformity';
inu.val       = {inu_reg,inu_co};
inu.help      = {'Specify the intensity nonuniformity (INU) settings for the current channel, which consist of a regularisation setting and a cutoff.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
label_files        = cfg_files;
label_files.tag    = 'images';
label_files.name   = 'Label maps';
label_files.filter = 'nifti';
label_files.num    = [1 Inf];
label_files.help   = {['Label maps are NIfTI images containing integer values, which must have the same ' ...
                       'dimensions and orientations as the scans of the corresponding subjects. ' ...
                       'Voxels of each value in the label map may be included in one or more tissue classes. ' ...
                       'For example, a label map showing the location of brain will include voxels that can ' ...
                       'be in grey or white matter classes. This information is specified in the ' ...
                       'confusion matrix.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cm_map        = cfg_entry;
cm_map.tag    = 'cm_map';
cm_map.name   = 'Row';
cm_map.strtype = 'n';
cm_map.num    = [1 Inf];
cm_map.help   = {'For this value in the label map, specify which tissue classes it can correspond to (including the K+1 implicit background class).'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cm        = cfg_repeat;
cm.tag    = 'cm';
cm.name   = 'Confusion matrix';
cm.values = {cm_map};
cm.help   = {'Specify rows of a confusion matrix, where each row corresponds to label values of 1, 2, ..., L + 1, etc in a label map.' ...
             'L are the number of labels in the label map. The last row (L + 1) needs to specify what classes unlabeled voxels can take.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%label_pr        = cfg_const;
%label_pr.tag    = 'w';
%label_pr.name   = 'Confidence';
%label_pr.val    = {0.99};
%label_pr.hidden = true;
%label_pr.help   = {'Degree of confidence in the labels.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
labels         = cfg_branch;
labels.tag     = 'true';
labels.name    = 'Has labels';
%labels.val    = {label_files,cm,const('w',0.99)};
labels.val     = {label_files};
labels.help    = {['If subjects have corresponding label maps to guide the ' ...
                   'segmentation, these need to be specified along with a ' ...
                   'confusion matrix that relates values in the label maps ' ...
                   'to which tissue classes they correspond with.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%no_labels      = cfg_const;
%no_labels.tag  = 'false';
%no_labels.name = 'No';
%no_labels.val  = {[]};
%no_labels.help = {'Subjects do not have corresponding label maps.'};
no_labels = const('false',[]);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
has_labels        = cfg_choice;
has_labels.tag    = 'labels';
has_labels.name   = 'Labels?';
has_labels.values = {labels,no_labels};
has_labels.val    = {no_labels};
has_labels.help   = {['Specify whether or not there are any pre-defined ' ...
                      'label maps for (all) the subjects in the current population.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
modality      = cfg_menu;
modality.tag  = 'modality';
modality.name = 'Modality';
modality.labels = {'MRI','CT'};
modality.values = {1,2};
modality.val    = {1};
modality.help   = {['Specify the modality of the scans in this channel. ' ...
                    'The main reason this is done is so that CT files can have ' ...
                    'a constant value of 1000 added to them to account for the ' ...
                    'way Hounsfield units are defined.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
chan          = cfg_branch;
chan.tag      = 'chan';
chan.name     = 'Channel';
chan.val      = {images, inu, modality};
chan.help     = {['There may be multiple scans of different modalities for each subject. ' ...
                  'These would be entered into different channels. Note that all scans ' ...
                  'within a channel should be the same modality.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
chans         = cfg_repeat;
chans.tag     = 'chans';
chans.name    = 'Channels';
chans.values  = {chan};
chans.num     = [1 Inf];
chans.val     = {chan};
chans.help    = {['Multiple image channels may be specified. For example, two channels may be used to contain ' ...
                  'the T2-weighted and PD-weighted scans of the subjects.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
pr_dat         = cfg_files;
pr_dat.tag     = 'file';
pr_dat.name    = 'Definition';
pr_dat.filter  = 'mat';
pr_dat.ufilter = '^prior.*\.mat$';
pr_dat.dir     = datadir;
pr_dat.num     = [0 1];
pr_dat.val     = {{}};
pr_dat.help    = {['Knowledge of Gaussian-Wishart priors for the intensity distributions of each cluster ' ...
                   'can help to inform the segmentation. When available, this information is specified in MATLAB prior*.m files. ' ...
                   'These files currently need to be hand-crafted. Unless you understand what you are doing, ' ...
                   'it is advised that you do not specify and intensity prior definition.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
pr_upd         = cfg_menu;
pr_upd.tag     = 'hyperpriors';
pr_upd.name    = 'Optimise';
pr_upd.labels  = {'Yes','No'};
pr_upd.values = {{'b0_priors',{0.01,0.01}}, []};
%pr_upd.values  = {{}, []};
pr_upd.val     = {pr_upd.values{2}};
pr_upd.help    = {['Specify whether the Gaussian-Wishart priors should be updated at each iteration. ' ...
                   'Enabling this can slow down convergence if there are small numbers of subjects. ' ...
                   'If only one subject is to be modelled (using a pre-computed template), then ' ...
                   'definitely turn off this option. However, if generating a template from a population of scans ' ...
                   'then this procedure is much more robust when the intensity priors are optimised.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
pr       = cfg_branch;
pr.tag   = 'pr';
pr.name  = 'Intensity prior';
pr.val   = {pr_dat,pr_upd};
pr.help  = {['Intensity distributions of each tissue class are modelled by a Gaussian distribution. ' ...
             'Prior knowledge about these distributions can make the model fitting more robust.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
pop       = cfg_branch;
pop.tag   = 'gmm';
pop.name  = 'Pop. of scans';
pop.val   = {chans, has_labels, pr,...
             const('tol_gmm', 0.0005), const('nit_gmm_miss',32), const('nit_gmm',8), const('nit_appear', 8)};
pop.check = @check_pop;
pop.help  = {'Information about a population of subjects that all have the same set of scans.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
pops            = cfg_repeat;
pops.tag        = 'pops';
pops.name       = 'Populations';
pops.values     = {pop};
pops.num        = [0 Inf];
pops.val        = {};
pops.help       = {['Multiple populations of subjects may be combined. For example, there may be ' ...
                    'T1-weighted scans and manually defined labels for one population, whereas ' ...
                    'another population may have T2-weighted and PD-weighted scans without labels. ' ...
                    'Yet another population might have CT scans. All subject''s data would be ' ...
                    'subdivided into the same tissue classes, although the intensity distributions ' ...
                    'of these tissues is likely to differ across populations.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
seg          = cfg_files;
seg.tag      = 'cat';
seg.name     = 'Class';
seg.filter   = 'nifti';
seg.ufilter  = '.*c[0-9].*';
seg.num      = [0 Inf];
seg.val      = {{}};
seg.help     = {['For each of the tissue class types, the user should specify the ' ...
                  'data to be included within the model fitting by selecting the files. ' ...
                  'It is important that the subject ordering of the files is the same ' ...
                  'across all classes.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
segs         = cfg_repeat;
segs.tag     = 'images';
segs.name    = 'Classes';
segs.values  = {seg};
segs.val     = {seg};
segs.help    = {['Images might have been segmented previously into a number of tissue classes. ' ...
                 'This framework allows such pre-segmented images to be included in the model fitting, ' ...
                 'in a similar way to the old Dartel toolbox for SPM. The user sets up a series ' ...
                 'of tissue class types (e.g., grey matter and white matter).']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
spop          = cfg_branch;
spop.tag      = 'cat';
spop.name     = 'Tissue class maps';
spop.val      = {images};
spop.check    = @check_segs;
spop.help     = {'UNUSED'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mu_exist         = cfg_files;
mu_exist.tag     = 'exist';
mu_exist.name    = 'Existing template';
mu_exist.filter  = 'nifti';
mu_exist.ufilter = '^mu_.*';
mu_exist.dir     = datadir;
mu_exist.num     = [1 1];
mu_exist.help    = {['The model can be fit using a previously computed template, which is not updated. ' ...
                     'Note that the template contains K-1 volumes within it, and that K should be compatible ' ...
                     'with various aspects of the data to which the model is fit. The template does not ' ...
                     'actually encode the tissue probabilities, but rather these probabilities can be generated ' ...
                     'from the template using a softmax function/* (${\bf p} = \frac{\exp {\boldsymbol{\mu}}}{1 + \sum_k \exp \mu_k}$)*/.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
nclass          = cfg_entry;
nclass.tag      = 'K';
nclass.name     = 'Number of classes';
nclass.strtype  = 'e';
nclass.val      = {9};
nclass.num      = [1 1];
nclass.help     = {['Specify K, the number of tissue classes encoded by the template. ' ...
                    'This value is ignored if it is incompatible with the specified data.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
vox          = cfg_entry;
vox.tag      = 'vx';
vox.name     = 'Voxel size';
vox.strtype  = 'e';
vox.val      = {1};
vox.num      = [1 1];
vox.help     = {['Specify the voxel size of the template to be created (mm). ' ...
                 'The algorithm will automatically attempt to determine suitable ' ...
                 'settings for its orientation and field of view.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mu_sett        = cfg_const;
mu_sett.tag    = 'mu_settings';
mu_sett.name   = 'Mu settings';
mu_sett.val    = {[1.0000e-05 0.5000 0]};
mu_sett.hidden = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mu_create       = cfg_branch;
mu_create.tag   = 'create';
mu_create.name  = 'Create template';
mu_create.val   = {nclass, vox, const('mu_settings',[1e-5 0.5 0]),const('issym',0)};
mu_create.help  = {['A tissue probability template will be constructed from all the aligned images. ' ...
                    'The algorithm alternates between re-computing the template and re-aligning all the ' ...
                    'images with this template. ' ...
                    'The user chooses the voxel size of the template and the number of tissue ' ...
                    'types it encodes.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mu_prov         = cfg_choice;
mu_prov.tag     = 'mu';
mu_prov.name    = 'Template';
mu_prov.values  = {mu_create, mu_exist};
mu_prov.val     = {mu_create};
mu_prov.help    = {[...
'The model can be run using a pre-computed template, or it can implicitly ' ...
'create an average shaped template from the population(s) of scan data. ' ...
'Here, the user chooses whether to create a template or use an existing one. ' ...
'Templates are named ``mu_*.nii``.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
aff             = cfg_menu;
aff.tag         = 'aff';
aff.name        = 'Affine';
aff.labels      = {'None', 'Translations only', 'Rigid body', 'Rigid+Zoom', 'Affine'};
aff.values      = {'', 'T(3)', 'SE(3)', 'SEZ(3)', 'Aff(3)'};
aff.val         = {'SE(3)'};
aff.help        = {[...
'Specify the type of affine transform to use in the model, which may be either ' ...
'none, translations only (T(3)) or rigid body (SE(3)). The fitting ' ...
'begins with affine registration, before continuing by interleaving ' ...
'affine and diffeomorphic registrations over multiple spatial scales.' ...
'Note that the "Affine" option is likely to throw up lots of ' ...
'warnings about ``QFORM0 representation has been rounded``, which ' ...
'can mostly be ignored.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dff             = cfg_entry;
dff.tag         = 'v_settings';
dff.name        = 'Shape regularisation';
dff.strtype     = 'e';
dff.num         = [1 5];
dff.val         = {[0.0001 0 0.4 0.1 0.4]};
dff.help        = {[...
'Specify the regularisation settings for the diffeomorphic registration. ' ...
'These consist of a vector of five values, which penalise different ' ...
'aspects of the warps:'],...
['    1. Absolute displacements need to be penalised by a tiny amount. ' ...
 'The first element encodes the amount of penalty on these. ' ...
 ' Ideally, absolute displacements should not be penalised, but it is usually ' ...
 ' necessary for technical reasons.'],...
['    2. The *membrane energy* of the deformation is penalised, usually by a ' ...
 'relatively small amount. This penalises the sum of squares of the ' ...
 'derivatives of the velocity field (i.e., the sum of squares of the elements ' ...
 'of the Jacobian tensors).'],...
['    3. The *bending energy* is penalised (3rd element). This penalises the ' ...
 'sum of squares of the 2nd derivatives of the velocity.'],...
['    4. Linear elasticity regularisation is also included. ' ...
 'This parameter (/*$\*/mu/*$*/) is similar to that for linear elasticity, ' ...
 'except it penalises the sum of squares of the Jacobian tensors after ' ...
 'they have been made symmetric (by averaging with the transpose). ' ...
 'This term essentially penalises length changes, without penalising rotations.'],...
['    5. The final term also relates to linear elasticity, and is the weight that ' ...
 'denotes how much to penalise changes to the divergence of the velocities (/*$\*/lambda/*$*/). ' ...
 'This divergence is a measure of the rate of volumetric expansion or contraction.'],...
'',...
'The default settings work reasonably well for most cases.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
odir            = cfg_files;
odir.tag        = 'odir';
odir.name       = 'Output directory';
odir.filter     = 'dir';
odir.num        = [1 1];
odir.val        = {{'.'}};
odir.help       = {'All output is written to the specified directory. If this is not specified, the current working directory will used by default.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
onam            = cfg_entry;
onam.tag        = 'onam';
onam.name       = 'Output name';
onam.strtype    = 's';
onam.val        = {'mb'};
onam.help       = {'Specify a key string for inclusion within all the output file names.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mb             = cfg_exbranch;
mb.tag         = 'run';
mb.name        = 'Fit Multi-Brain model';
mb.val         = {mu_prov, aff, dff, const('del_settings',Inf), onam, odir, segs, pops,...
                   const('accel',0.8), const('min_dim', 8), const('tol',0.001),...
                   const('sampdens',2),const('save',true),const('nworker',0)};
mb.prog        = @run_mb;
mb.vout        = @vout_mb_run;
mb.check       = @check_images;
mb.help        = {...
    ['This is where the model fitting is actually done. ' ...
     'The outputs from the model fitting are an initial velocity field ' ...
     'and deformation field for each subject. In addition, fitting the ' ...
     'model may also generate a representation of a set of average shaped ' ...
     'tissue probability maps that can serve as a template for future ' ...
     'model fitting. ',...
     'The model /* \cite{brudfors2020flexible} */ builds on a number of previous works and has ' ...
     'the objective of achieving diffeomorphic alignment of a wide variety of medical ' ...
     'image modalities into a common anatomical space. ' ...
     'This involves the ability to construct a "tissue probability template" from ' ...
     'one or more populations of scans through group-wise alignment/* \cite{john_averageshape,blaiotta2018generative}*/. '], ...
    ['Diffeomorphic deformations are computed within a geodesic shooting framework' ...
     '/* \cite{ashburner2011diffeomorphic}*/, which is optimised with a Gauss-Newton ' ...
     'strategy that uses a multi-grid approach to solve the system of linear equations/* \cite{ashburner07}*/. ' ...
     'Variability among image contrasts is modelled using a more sophisticated version of ' ...
     'the Gaussian mixture model with bias correction framework originally proposed ' ...
     'in the "Unified Segmentation" paper/* \cite{ashburner05}*/, ' ...
     'and which has been extended to account for known variability of the intensity ' ...
     'distributions of different tissues/* \cite{blaiotta2016variational}*/. ' ...
     'This model has been shown to provide a good model of the intensity distributions ' ...
     'of different imaging modalities/* \cite{brudfors2019empirical}*/.'],...
     'This work was funded by the EU Human Brain Project''s Grant Agreement No 785907 (SGA2).'};
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
res_file         = cfg_files;
res_file.tag     = 'result';
res_file.name    = 'MB results file';
res_file.filter  = 'mat';
res_file.ufilter = '^mb_fit.*';
res_file.num     = [1 1];
res_file.help    = {['Specify the results file obtained from running **Fit Multi-Brain model**. ' ...
                     'This will be named ``mb_fit_*.mat`` and contain a link to where any resulting ' ...
                     'template may be found.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
ix        = cfg_entry;
ix.tag    = 'ix';
ix.name   = 'Indices';
ix.strtype = 'n';
ix.num    = [1 Inf];
ix.help   = {['Specify indices. For example, if the original model ' ...
              'had K=9 and you wish to combine the final three classes, ' ...
              'then enter ``1 2 3 4 5 6 7 7 7``. Note that K refers to the ' ...
              'total number of tissue maps -- including the implicit ' ...
              'background.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {['The bounding box (in voxels) of the merged template volume. ' ...
               'Non-finite values indicate to use original template dimensions.']};
bb.strtype = 'r';
bb.num     = [2 3];
bb.val     = {[1 1 1; Inf Inf Inf]};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Uses previously defined onam
onam.val = {'merged'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mrg      = cfg_exbranch;
mrg.tag  = 'merge';
mrg.name = 'Merge tissues';
mrg.val  = {res_file, ix, bb, onam, odir};
mrg.prog = @spm_mb_merge;
mrg.vout = @vout_mb_merge;
mrg.help = {['This option is for merging template tissues together and ' ...
             'extracting intensity priors for later use when running ' ...
             '**Fit Multi-Brain model**. This is typically used for ' ...
             're-ordering tissue classes or combining multiple classes ' ...
             '(e.g. from air, which has a non-Gaussian intensity distribution ' ...
             'that is often has multiple "tissues" fitted to it) ' ...
             'into one. Generated templates also usually have ' ...
             'a large field of view, so it is often desirable to trim them down ' ...
             'so the field of view covers a smaller region of anatomy.']};
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
res_file         = cfg_files;
res_file.tag     = 'result';
res_file.name    = 'MB results file';
res_file.filter  = 'mat';
res_file.ufilter = '^mb_fit.*';
res_file.num     = [1 1];
res_file.help    = {'Specify the results file obtained from previously running **Fit Multi-Brain model**.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
i        = cfg_menu;
i.tag    = 'i';
i.name   = 'Images';
i.labels = {'No','Yes'};
i.values = {false,true};
i.val    = {false};
i.help   = {['Specify whether versions of the original images, ' ...
             'but with missing values filled in, should be written out.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mi        = cfg_menu;
mi.tag    = 'mi';
mi.name   = 'INU corrected';
mi.labels = {'No','Yes'};
mi.values = {false,true};
mi.val    = {false};
mi.help   = {['Specify whether INU corrected versions of the original images ' ...
              '(with missing values filled in) should be written out.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
wi        = cfg_menu;
wi.tag    = 'wi';
wi.name   = 'Warped images';
wi.labels = {'No','Yes'};
wi.values = {false,true};
wi.val    = {false};
wi.help   = {['Specify whether spatially normalised versions of the images ' ...
              '(missing values filled in) should be written out.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
wmi        = cfg_menu;
wmi.tag    = 'wmi';
wmi.name   = 'Warped INU corrected';
wmi.labels = {'No','Yes'};
wmi.values = {false,true};
wmi.val    = {false};
wi.help   = {['Specify whether spatially normalised versions of the INU corrected images ' ...
              '(missing values filled in) should be written out.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
inu        = cfg_menu;
inu.tag    = 'inu';
inu.name   = 'INU';
inu.labels = {'No', 'Yes'};
inu.values = {false,true};
inu.val    = {false};
inu.help   = {'Specify whether the estimated INU fields should be written out.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Tissues';
c.strtype = 'n';
c.num     = [0 Inf];
c.val     = {[]};
c.help    = {'Specify the indices of any native-space tissue class images to be written.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
wc         = cfg_entry;
wc.tag     = 'wc';
wc.name    = 'Warped tissues';
wc.strtype = 'n';
wc.num     = [0 Inf];
wc.val     = {[]};
wc.help    = {'Specify the indices of any spatially normalised tissue class images to be written.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mwc         = cfg_entry;
mwc.tag     = 'mwc';
mwc.name    = 'Warped mod. tissues';
mwc.strtype = 'n';
mwc.num     = [0 Inf];
mwc.val     = {[]};
mwc.help    = {'Specify the indices of any spatially normalised and Jacobian-scaled ("modulated") tissue class images to be written.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
sm         = cfg_entry;
sm.tag     = 'sm';
sm.name    = 'Scalar momentum';
sm.strtype = 'n';
sm.num     = [0 Inf];
sm.val     = {[]};
sm.help    = {'Specify the indices of any scalar momentums to be written.'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mrf         = cfg_entry;
mrf.tag     = 'mrf';
mrf.name    = 'MRF Parameter';
mrf.help    = {'When tissue class images are written out, a few iterations of a simple Markov random field (MRF) cleanup procedure are run.  This parameter controls the strength of the MRF. Setting the value to zero will disable the cleanup.'};
mrf.strtype = 'r';
mrf.num     = [1 1];
mrf.val     = {0};
mrf.hidden  = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Output smoothing';
fwhm.strtype = 'e';
fwhm.num     = [1 3];
fwhm.val     = {[0 0 0]};
fwhm.help    = {['Full width at half maximum (FWHM) of Gaussian smoothing kernel ' ...
                 'for smoothing of template space data: warped tissues, modulated warped tissues ' ...
                 'and scalar momentums.']};
fwhm.hidden  = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {['The bounding box of the template space output, as a 2x3 array. ' ...
               'The specified values should be the minimum and maximum x coordinates, ' ...
               'the minimum and maximum y coordinates, and the minimum and maximum z coordinates.']};
bb.strtype = 'r';
bb.num     = [2 3];
bb.val     = {[NaN NaN NaN; NaN NaN NaN]};
bb.hidden  = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel sizes';
vox.help    = {'The (isotropic) voxel size of the template space output, in mm.'};
vox.strtype = 'r';
vox.num     = [1 1];
vox.val     = {NaN};
vox.hidden  = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
proc_zn        = cfg_entry;
proc_zn.tag    = 'proc_zn';
proc_zn.name   = 'Process responsibilities';
proc_zn.help   = {['Function for processing native space responsibilities, ' ...
                  'given as a function handle ``@(x) foo(x)``. The argument (``x``) is of ' ...
                  '``size(x) = [1, 4]``, where the first three dimensions are the size ' ...
                  'of the image and the last dimension is the number of segmentation ' ...
                  'classes (K + 1).']};
proc_zn.val    = {{}};
proc_zn.hidden = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
odir        = cfg_files;
odir.tag    = 'odir';
odir.name   = 'Output directory';
odir.filter = 'dir';
odir.num    = [1 1];
odir.val    = {{''}};
odir.help   = {'All output is written to the specified directory. If this is not specified, the output directory of the run module will be used by default.'};
odir.hidden = true;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
out      = cfg_exbranch;
out.tag  = 'out';
out.name = 'Output';
out.val  = {res_file, i, mi, wi, wmi, inu, c, wc, mwc, sm, mrf, fwhm, bb, vox, proc_zn, odir};
out.prog = @spm_mb_output;
out.help = {[...
'When **Fit Multi-Brain model** is run, the resulting model fit contains ' ...
'information that allows a lot of other derived images to be generated. ' ...
'For example, the results file encodes the INU fields, which allows ' ...
'images to be INU corrected. ' ...
'It contains information about intensity distributions, which (when ' ...
'combined with other extracted information) enables tissue segmentation ' ...
'to be achieved. ' ...
'And of course, the estimated deformations allow spatially normalised versions ' ...
'of these results to be generated.'], [...
'The **Output** functionality allows a range of derived data to be ' ...
'generated from these parameter estimates. When it is executed, it generates ' ...
'the derived data from the images originally entered into **Fit Multi-Brain model**. ' ...
'In order to do this, it needs to assume that the images, as well as the ' ...
'results files, have not been moved from their original locations. ' ...
'If they can not be fund, then **Output** will crash out in a not very elegant way.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
images        = cfg_files;
images.tag    = 'images';
images.name   = 'Scans';
images.filter = 'nifti';
images.num    = [1 Inf];
images.help   = {['Select one NIfTI format scan for each subject.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
fil      = cfg_exbranch;
fil.tag  = 'fil';
fil.name = 'Image Labelling';
fil.val  = {images};
fil.prog = @spm_label;
fil.vout = @vout_fil;
fil.help = {...
    'Factorisation-based Image Labelling.', ...
    ['Label brains according to the Neuromorphometrics protocol using the method '...
     'of Yan et al/* \cite{yan2021factorisation} */. Note that the algorithm '...
     'may need to download additional data from Figshare:'], ...
     '    * https://figshare.com/articles/dataset/Trained_FIL_model/17143370/2',...
     '    * https://figshare.com/articles/dataset/Head_Tissue_Template/17143289'};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
images        = cfg_files;
images.tag    = 'images';
images.name   = 'Scans';
images.filter = 'nifti';
images.num    = [1 Inf];
images.help   = {['Select one NIfTI format scan for each subject.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
mbnorm      = cfg_exbranch;
mbnorm.tag  = 'mbnorm';
mbnorm.name = 'Spatially normalise using Multi-Brain';
mbnorm.val  = {images};
mbnorm.prog = @spm_mbnorm;
mbnorm.vout = @vout_mbnorm;
mbnorm.help = {...
    ['Spatially normalise by matching a pre-computed template to each image. '...
    'The template is not in ICBM space, but a deformation is provided for composing '...
    'with the inverse of the estimated deformations in order to obtain a suitable '...
    'spatially normalising deformation. '...
    'Note that the algorithm '...
    'may need to download additional data from Figshare:'], '', ...
    '    * https://figshare.com/articles/dataset/Deformation_field_mapping_from_ICMB_to_mu_X_nii/17143796/1',...
    '    * https://figshare.com/articles/dataset/Head_Tissue_Template/17143289', ...
    '', ...
    'This approach is described in:', '', ...
    ['    * Brudfors, M., Balbastre, Y., Flandin, G., Nachev, P. and Ashburner, J., 2020. ' ...
     'Flexible Bayesian modelling for nonlinear image registration. ' ...
     'In Medical Image Computing and Computer Assisted Intervention–MICCAI 2020: ' ...
     '23rd International Conference, Lima, Peru, October 4–8, 2020, ' ...
     'Proceedings, Part III 23 (pp. 253-263). Springer International Publishing.'],...
    ['    * Blaiotta, C., Freund, P., Cardoso, M.J. and Ashburner, J., 2018. ' ...
     'Generative diffeomorphic modelling of large MRI data sets for probabilistic template construction. ' ...
     'NeuroImage, 166, pp.117-134.']};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cfg        = cfg_choice;
cfg.tag    = 'mb';
cfg.name   = 'Multi-Brain toolbox';
cfg.values = {mb,mrg,out,fil,mbnorm};
cfg.help   = {['The Multi-Brain (MB) toolbox has the general aim of integrating a number of disparate ' ...
               'image analysis components within a single unified generative modelling framework ' ...
               '(segmentation, nonlinear registration, image translation, etc.).']};
%__________________________________________________________________________
%
%__________________________________________________________________________
function cfg = const(tag,val)
cfg          = cfg_const;
cfg.tag      = tag;
cfg.val      = {val};
cfg.hidden   = true;
%__________________________________________________________________________
%
%__________________________________________________________________________
function  out = run_mb(cfg)
[dat,sett]    = spm_mb_init(cfg);
if ~isempty(dat)
    [dat,sett] = spm_mb_fit(dat,sett);
    out        = out_mb_run(sett,dat);
    save(out.fit{1},'sett','dat');
else
    out           = struct('fit',{{''}},'mu',{{''}},'v',{{''}},'psi',{{''}});
    warning('No images specified: Done nothing.');
end
%__________________________________________________________________________
%
%__________________________________________________________________________
function dr = datadir
fullname    = mfilename('fullpath');
pth         = fileparts(fullname);
dr          = fullfile(pth,'data');
%__________________________________________________________________________
%
%__________________________________________________________________________
function str = check_images(cfg)
has_data = false;
for c=1:numel(cfg.gmm)
    if numel(cfg.gmm(c).chan)>0
        if numel(cfg.gmm(c).chan(1).images)>0
            has_data = true;
            break
        end
    end
end
if numel(cfg.cat)>0
    if numel(cfg.cat{1})>0
        if ~isempty(cfg.cat{1}{1})
            has_data = true;
        end
    end
end
if has_data
    str = {};
else
    str = {'No images specified.'};
end
%__________________________________________________________________________
%
%__________________________________________________________________________
function str = check_pop(cfg)
str = {};
N   = -1;
for c=1:numel(cfg.chan)
    Nc = numel(cfg.chan(c).images);
    if N<0
        N = Nc;
    else
        if N~=Nc
            str = {'Incompatible numbers of scans over channels.'};
        end
    end
end
if isfield(cfg.labels,'true')
    Nc = numel(cfg.chan(c).images);
    if N~=Nc
        str = {'Incompatible numbers of label images.'};
    end
end
%__________________________________________________________________________
%
%__________________________________________________________________________
function str = check_segs(cfg)
str    = {};
N      = -1;
images = cfg.images;
for c=1:numel(images)
    Nc = numel(images{c});
    if N<0
        N = Nc;
    else
        if N~=Nc
            str = {'Incompatible numbers of categorical images.'};
        end
    end
end
%__________________________________________________________________________
%
%__________________________________________________________________________
function dep = vout_mb_run(cfg)

if isfield(cfg.mu,'create')
    mudep            = cfg_dep;
    mudep.sname      = 'Template';
    mudep.src_output = substruct('.','mu','()',{':'});
    mudep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
else
    mudep = {};
end

matdep            = cfg_dep;
matdep.sname      = 'MB-fit';
matdep.src_output = substruct('.','fit','()',{':'});
matdep.tgt_spec   = cfg_findspec({{'filter','mat'}});

vdep            = cfg_dep;
vdep.sname      = 'Velocity Fields';
vdep.src_output = substruct('.','v','()',{':'});
vdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

ydep            = cfg_dep;
ydep.sname      = 'Deformation Fields';
ydep.src_output = substruct('.','psi','()',{':'});
ydep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

dep = [matdep mudep vdep ydep];
%__________________________________________________________________________
%
%__________________________________________________________________________
function out = out_mb_run(sett,dat)
if isfield(sett.mu,'create')
    out = struct('fit',[],'mu',[],'v',[],'psi',[]);
    if isa(sett.mu.create.mu,'nifti')
        out.mu = {sett.mu.create.mu.dat.fname};
    else
        out.mu = {sett.mu.create.mu};
    end
else
    out = struct('fit',[],'v',[],'psi',[]);
end
out.fit = {fullfile(sett.odir,['mb_fit_' sett.onam '.mat'])};
out.v   = cell(numel(dat),1);
out.psi = cell(numel(dat),1);
for n=1:numel(dat)
    out.v{n}   = dat(n).v.dat.fname;
    out.psi{n} = dat(n).psi.dat.fname;
end
%__________________________________________________________________________
%
%__________________________________________________________________________
function dep = vout_mb_merge(cfg)
matdep            = cfg_dep;
matdep.sname      = 'Intensity Priors';
matdep.src_output = substruct('.','priors','()',{':'});
matdep.tgt_spec   = cfg_findspec({{'filter','mat'}});

mudep            = cfg_dep;
mudep.sname      = 'Merged Template';
mudep.src_output = substruct('.','mu','()',{':'});
mudep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

dep = [mudep, matdep];
%__________________________________________________________________________
%
%__________________________________________________________________________
function dep = vout_fil(cfg)
dep            = cfg_dep;
dep.sname      = 'Labelled brains';
dep.src_output = substruct('.','labels','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%__________________________________________________________________________
%
%__________________________________________________________________________
function dep = vout_mbnorm(cfg)
vo1            = cfg_dep;
vo1.sname      = 'Deformation';
vo1.src_output = substruct('.','def');
vo1.tgt_spec   = cfg_findspec({{'filter','nifti'}});
vo2            = cfg_dep;
vo2.sname      = 'Warped Images';
vo2.src_output = substruct('.','warped');
vo2.tgt_spec   = cfg_findspec({{'filter','image'}});
dep            = [vo1, vo2];

