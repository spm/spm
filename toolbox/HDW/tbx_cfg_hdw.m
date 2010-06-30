function hdw = tbx_cfg_hdw
% Configuration file for toolbox 'High-Dimensional Warping'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_hdw.m 3960 2010-06-30 17:41:24Z ged $

% ---------------------------------------------------------------------
% ref Reference Image
% ---------------------------------------------------------------------
ref         = cfg_files;
ref.tag     = 'ref';
ref.name    = 'Reference Image';
ref.help    = {'This is the reference image, which remains stationary.'};
ref.filter  = 'image';
ref.ufilter = '.*';
ref.num     = [1 1];
% ---------------------------------------------------------------------
% mov Moved Image
% ---------------------------------------------------------------------
mov         = cfg_files;
mov.tag     = 'mov';
mov.name    = 'Moved Image';
mov.help    = {'This is the moved image, which is warped to match the reference.'};
mov.filter  = 'image';
mov.ufilter = '.*';
mov.num     = [1 1];
% ---------------------------------------------------------------------
% data Subject
% ---------------------------------------------------------------------
data1       = cfg_branch;
data1.tag   = 'data';
data1.name  = 'Subject';
data1.val   = {ref mov};
data1.help  = {'Two images of the same subject, which are to be registered together.  Prior to nonlinear high-dimensional warping, the images should be rigidly registered with each other.'};
% ---------------------------------------------------------------------
% data Subjects
% ---------------------------------------------------------------------
data         = cfg_repeat;
data.tag     = 'data';
data.name    = 'Subjects';
data.help    = {'Specify pairs of images to match together.'};
data.values  = {data1};
data.num     = [0 Inf];
% ---------------------------------------------------------------------
% nits Iterations
% ---------------------------------------------------------------------
nits         = cfg_entry;
nits.tag     = 'nits';
nits.name    = 'Iterations';
nits.help    = {'Number of iterations for the bias correction'};
nits.val     = {8};
nits.strtype = 'n';
nits.num     = [1  1];
% ---------------------------------------------------------------------
% fwhm Bias FWHM
% ---------------------------------------------------------------------
fwhm         = cfg_menu;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Bias FWHM';
fwhm.help    = {'FWHM of Gaussian smoothness of bias. If your intensity nonuniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity nonuniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity nonuniformities.'};
fwhm.val     = {60};
fwhm.labels = {'30mm cutoff'
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
               'No correction'}';
fwhm.values = {30 40 50 60 70 80 90 100 110 120 130 140 150 Inf};
% ---------------------------------------------------------------------
% reg Bias regularisation
% ---------------------------------------------------------------------
reg         = cfg_menu;
reg.tag     = 'reg';
reg.name    = 'Bias regularisation';
reg.help    = {'We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity nonuniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of nonuniformity.'};
reg.val     = {1e-6};
reg.labels = {'no regularisation'
              'extremely light regularisation'
              'very light regularisation'
              'light regularisation'
              'medium regularisation'
              'heavy regularisation'
              'very heavy regularisation'
              'extremely heavy regularisation'}';
reg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
% ---------------------------------------------------------------------
% lmreg Levenberg-Marquardt regularisation
% ---------------------------------------------------------------------
lmreg         = cfg_menu;
lmreg.tag     = 'lmreg';
lmreg.name    = 'Levenberg-Marquardt regularisation';
lmreg.help    = {'Levenberg-Marquardt regularisation keeps the bias correction part stable. Higher values means more stability, but slower convergence.'};
lmreg.val     = {1e-6};
lmreg.labels  ={'no regularisation'
                'extremely light regularisation'
                'very light regularisation'
                'light regularisation'
                'medium regularisation'
                'heavy regularisation'
                'very heavy regularisation'
                'extremely heavy regularisation'
                }';
lmreg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
% ---------------------------------------------------------------------
% bias_opts Bias Correction Options
% ---------------------------------------------------------------------
bias_opts       = cfg_branch;
bias_opts.tag   = 'bias_opts';
bias_opts.name  = 'Bias Correction Options';
bias_opts.val   = {nits fwhm reg lmreg };
bias_opts.help  = {'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   'Before registering the images, an approximate bias correction is estimated for the moved image. This is based on minimising the difference between the images an a symmetric way. Prior to registering the images, they should be rigidly aligned together.  The bias correction is estimated once for these aligned images.'}';
% ---------------------------------------------------------------------
% nits Iterations
% ---------------------------------------------------------------------
nitsw         = cfg_entry;
nitsw.tag     = 'nits';
nitsw.name    = 'Iterations';
nitsw.help    = {'Number of iterations for the warping.'};
nitsw.val     = {8};
nitsw.strtype = 'n';
nitsw.num     = [1  1];
% ---------------------------------------------------------------------
% reg Warping regularisation
% ---------------------------------------------------------------------
regw         = cfg_entry;
regw.tag     = 'reg';
regw.name    = 'Warping regularisation';
regw.help    = {'There is a tradeoff between the smoothness of the estimated warps, and the difference between the registered images.  Higher values mean smoother warps, at the expense of a lower mean squared difference between the images.'};
regw.val     = {4};
regw.strtype = 'e';
regw.num     = [1  1];
% ---------------------------------------------------------------------
% warp_opts Warping Options
% ---------------------------------------------------------------------
warp_opts      = cfg_branch;
warp_opts.tag  = 'warp_opts';
warp_opts.name = 'Warping Options';
warp_opts.val  = {nitsw regw};
warp_opts.help = {'There are a couple of user-customisable warping options.'};
% ---------------------------------------------------------------------
% hdw High-Dimensional Warping
% ---------------------------------------------------------------------
hdw         = cfg_exbranch;
hdw.tag     = 'hdw';
hdw.name    = 'High-Dimensional Warping';
hdw.val     = {data bias_opts warp_opts};
hdw.help    = {'This toolbox is a Bayesian method for three dimensional registration of brain images/* \cite{ashburner00a} */. A finite element approach is used to obtain a maximum a posteriori (MAP) estimate of the deformation field at every voxel of a template volume.  The priors used by the MAP estimate penalize unlikely deformations and enforce a continuous one-to-one mapping.  The deformations are assumed to have some form of symmetry, in that priors describing the probability distribution of the deformations should be identical to those for the inverses (i.e., warping brain A to brain B should not be different probablistically from warping B to A).  A gradient descent algorithm is used to estimate the optimum deformations.'
               'Deformation fields are written with the same name as the moved image, but with "y_" prefixed on to the filename.  Jacobian determinant images are also written (prefixed by "jy_").'}';
hdw.prog    = @spm_local_hdw;

%======================================================================
function spm_local_hdw(job)

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','HDW')); end
spm_hdw(job);