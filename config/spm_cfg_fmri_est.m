function fmri_est = spm_cfg_fmri_est
% SPM Configuration file for Model Estimation
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_fmri_est.m 5959 2014-04-16 17:14:33Z will $


%==========================================================================
% spmmat Select SPM.mat
%==========================================================================
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {
                  'Select the SPM.mat file that contains the design specification. '
                  'The directory containing this file is known as the input directory.'
}';
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

%==========================================================================
% Classical Classical Inference
%==========================================================================
Classical      = cfg_const;
Classical.tag  = 'Classical';
Classical.name = 'Classical';
Classical.val  = {1};
Classical.help = {
                     'Model parameters are estimated using Restricted Maximum Likelihood (ReML). This assumes the error correlation structure is the same at each voxel. This correlation can be specified using either an AR(1) or an Independent and Identically Distributed (IID) error model. These options are chosen at the model specification stage. ReML estimation should be applied to spatially smoothed functional images.'
                     ''
                     'After estimation, specific profiles of parameters are tested using a linear compound or contrast with the T or F statistic. The resulting statistical map constitutes an SPM. The SPM{T}/{F} is then characterised in terms of focal or regional differences by assuming that (under the null hypothesis) the components of the SPM (ie. residual fields) behave as smooth stationary Gaussian fields.'
}';

%--------------------------------------------------------------------------
% volBlocktype Block type
%--------------------------------------------------------------------------
volBlocktype        = cfg_menu;
volBlocktype.tag    = 'block_type';
volBlocktype.name   = 'Block type';
volBlocktype.val    = {'Slices'};
volBlocktype.help   = {'Enter the block type, i.e. "Slices" or "Subvolumes"'}';
volBlocktype.labels = {'Slices', 'Subvolumes'}';
volBlocktype.values = {'Slices', 'Subvolumes'}';

%--------------------------------------------------------------------------
% Volume Volume
%--------------------------------------------------------------------------
Volume      = cfg_branch;
Volume.tag  = 'volume';
Volume.name = 'Volume';
Volume.val  = {volBlocktype};
Volume.help = {'A volume of data is analysed in "blocks", which can be a slice or 3D subvolume, where the extent of each subvolume is determined using a graph partitioning algorithm. Enter the block type, i.e. "Slices" or "Subvolumes".'};

%--------------------------------------------------------------------------
% Slices Slices
%--------------------------------------------------------------------------
SliceNs         = cfg_entry;
SliceNs.tag     = 'numbers';
SliceNs.name    = 'Slice numbers';
SliceNs.help    = {' '};
SliceNs.strtype = 'n';
SliceNs.num     = [Inf 1];

%--------------------------------------------------------------------------
% slBlocktype Block type
%--------------------------------------------------------------------------
slBlocktype        = cfg_menu;
slBlocktype.tag    = 'block_type';
slBlocktype.name   = 'Block type';
slBlocktype.val    = {'Slices'};
slBlocktype.help   = {'Enter the block type, i.e. "Slices" or "Subvolumes"'}';
slBlocktype.labels = {'Slices', 'Subvolumes'}';
slBlocktype.values = {'Slices', 'Subvolumes'}';

%--------------------------------------------------------------------------
% Slices Slices
%--------------------------------------------------------------------------
Slices      = cfg_branch;
Slices.tag  = 'slices';
Slices.name = 'Slices';
Slices.val  = {SliceNs slBlocktype};
Slices.help = {'Enter Slice Numbers. This can be a single slice or multiple slices. If you select a single slice or only a few slices you must be aware of the interpolation options when, after estimation, displaying the estimated images eg. images of contrasts or AR maps. The default interpolation option may need to be changed to nearest neighbour (NN) (see bottom right hand of graphics window) for you slice maps to be visible.'};

%--------------------------------------------------------------------------
% Clustermask Cluster mask
%--------------------------------------------------------------------------
Clustermask         = cfg_files;
Clustermask.tag     = 'mask';
Clustermask.name    = 'Cluster mask';
Clustermask.help    = {'Select cluster image'}';
Clustermask.filter  = 'image';
Clustermask.ufilter = '.*';
Clustermask.num     = [0 1];

%--------------------------------------------------------------------------
%  clBlocktype Block type
%--------------------------------------------------------------------------
clBlocktype        = cfg_menu;
clBlocktype.tag    = 'block_type';
clBlocktype.name   = 'Block type';
clBlocktype.val    = {'Slices'};
clBlocktype.help   = {'Enter the block type, i.e. "Slices" or "Subvolumes"'}';
clBlocktype.labels = {'Slices', 'Subvolumes'}';
clBlocktype.values = {'Slices', 'Subvolumes'}';

%--------------------------------------------------------------------------
% Clusters Clusters
%--------------------------------------------------------------------------
Clusters      = cfg_branch;
Clusters.tag  = 'clusters';
Clusters.name = 'Clusters';
Clusters.val  = {Clustermask clBlocktype};
Clusters.help = {'Because estimation can be time consuming an option is provided to analyse selected clusters rather than the whole volume.'};

%--------------------------------------------------------------------------
% space Analysis Space
%--------------------------------------------------------------------------
space         = cfg_choice;
space.tag     = 'space';
space.name    = 'Analysis Space';
space.val     = {Volume};
space.help    = {'Because estimation can be time consuming options are provided to analyse selected slices or clusters rather than the whole volume.'};
space.values  = {Volume Slices Clusters};

%--------------------------------------------------------------------------
% LogEv - Compute F 
%--------------------------------------------------------------------------
LogEv         = cfg_menu;
LogEv.tag     = 'LogEv';
LogEv.name    = 'Log evidence map';
LogEv.val     = {'No'};
LogEv.help    = {
                'Computes the log evidence for each voxel. When comparing models that differ, for example, just in their'
                'parametric modulation columns, it is important that these columns should be identically scaled.'
                'That is, they should have the same variance. If you just wish to compare models over a local region'
                'it may be better to use the command line option spm_vb_regionF.m. This latter function automatically scales'
                'design matrix columns (except for the mean) to have unit variance'};
LogEv.labels  = {'No','Yes'}';
LogEv.values  = {'No','Yes'}';

%--------------------------------------------------------------------------
% signal Signal priors
%--------------------------------------------------------------------------
signal         = cfg_menu;
signal.tag     = 'signal';
signal.name    = 'Signal priors';
signal.help    = {
                  '[UGL] Unweighted Graph Laplacian. This spatial prior is the recommended option. Regression coefficients at a given voxel are (softly) constrained to be similar to those at nearby voxels. The strength of this constraint is determined by a spatial precision parameter that is estimated from the data. Different regression coefficients have different spatial precisions allowing each putative experimental effect to have its own spatial regularity. '
                  ''
                  '[GMRF] Gaussian Markov Random Field. This is equivalent to a normalized UGL. '
                  ''
                  '[LORETA] Low resolution Tomography Prior. This is equivalent to UGL squared. It is a standatd choice for EEG source localisation algorithms. '
                  ''
                  '[WGL] Weighted Graph Laplacian. This is a generalization of the UGL, where weights can be used to preserve "edges" of functional responses.'
                  ''
                  '[Global] Global Shrinkage prior. This is not a spatial prior in the sense that regression coefficients are constrained to be similar to neighboring voxels. Instead, the average effect over all voxels (global effect) is assumed to be zero and all regression coefficients are shrunk towards this value in proporation to the prior precision. This is the same prior that is used for Bayesian estimation at the second level models, except that here the prior precision is estimated separaetly for each slice. '
                  ''
                  '[Uninformative] A flat prior. Essentially, no prior information is used. If you select this option then VB reduces to Maximum Likelihood (ML)estimation. This option is useful if, for example, you do not wish to use a spatial prior but wish to take advantage of the voxel-wise AR(P) modelling of noise processes. In this case, you would apply the algorithm to images that have been spatially smoothed. For P=0, ML estimation in turn reduces to Ordinary Least Squares (OLS) estimates, and for P>0 ML estimation is equivalent to a weighted least squares (WLS) but where the weights are different at each voxel (reflecting the different noise correlation at each voxel). '
}';
signal.labels = {
                 'UGL'
                 'GMRF'
                 'LORETA'
                 'WGL'
                 'Global'
                 'Uninformative'
}';
signal.values = {
                 'UGL'
                 'GMRF'
                 'LORETA'
                 'WGL'
                 'Global'
                 'Uninformative'
}';
signal.val     = {'UGL'};

%--------------------------------------------------------------------------
% ARP AR model order
%--------------------------------------------------------------------------
ARP         = cfg_entry;
ARP.tag     = 'ARP';
ARP.name    = 'AR model order';
ARP.help    = {
               'An AR model order of 3 is the default. Cardiac and respiratory artifacts are periodic in nature and therefore require an AR order of at least 2. In previous work, voxel-wise selection of the optimal model order showed that a value of 3 was the highest order required. '
               ''
               'Higher model orders have little effect on the estimation time. If you select a model order of zero this corresponds to the assumption that the errors are IID. This AR specification overrides any choices that were made in the model specification stage.'
               ''
               'Voxel-wise AR models are fitted separately for each session of data. For each session this therefore produces maps of AR(1), AR(2) etc coefficients in the output directory. '
}';
ARP.strtype = 'n';
ARP.num     = [Inf 1];
ARP.val     = {3};

%--------------------------------------------------------------------------
% UGL UGL
%--------------------------------------------------------------------------
UGL         = cfg_const;
UGL.tag     = 'UGL';
UGL.name    = 'UGL';
UGL.val     = {1};
UGL.help    = {'[UGL] Unweighted graph-Laplacian. This is the default option. This spatial prior is the same as that used for the regression coefficients. Spatial precisions are estimated separately for each AR coefficient eg. the AR(1) coefficient over space, AR(2) over space etc. '};

%--------------------------------------------------------------------------
% GMRF GMRF
%--------------------------------------------------------------------------
GMRF        = cfg_const;
GMRF.tag    = 'GMRF';
GMRF.name   = 'GMRF';
GMRF.val    = {1};
GMRF.help   = {'[GMRF] Gaussian Markov Random Field. See comments on GMRF priors for regresion coefficients. '};

%--------------------------------------------------------------------------
% LORETA LORETA
%--------------------------------------------------------------------------
LORETA      = cfg_const;
LORETA.tag  = 'LORETA';
LORETA.name = 'LORETA';
LORETA.val  = {1};
LORETA.help = {'[LORETA] Low resolution Tomography Prior. See comments on LORETA priors for regresion coefficients.'};

%--------------------------------------------------------------------------
% tissue_type Tissue-type
%--------------------------------------------------------------------------
tissue_type         = cfg_files;
tissue_type.tag     = 'tissue_type';
tissue_type.name    = 'Tissue-type';
tissue_type.help    = {
                       '[Tissue-type] AR estimates at each voxel are biased towards typical values for that tissue type (eg. gray, white, CSF). If you select this option you will need to then select files that contain tissue type maps (see below). These are typically chosen to be Grey Matter, White Matter and CSF images derived from segmentation of registered structural scans.'
                       ''
                       'Previous work has shown that there is significant variation in AR values with tissue type. However, GMRF priors have previously been favoured by Bayesian model comparison.'
}';
tissue_type.filter  = 'image';
tissue_type.ufilter = '.*';
tissue_type.num     = [1 Inf];

%--------------------------------------------------------------------------
% Robust Robust
%--------------------------------------------------------------------------
Robust         = cfg_const;
Robust.tag     = 'Robust';
Robust.name    = 'Robust';
Robust.val     = {1};
Robust.help    = {'Robust GLM. Uses Mixture of Gaussians noise model.'};

%--------------------------------------------------------------------------
% noise Noise priors
%--------------------------------------------------------------------------
noise         = cfg_choice;
noise.tag     = 'noise';
noise.name    = 'Noise priors';
noise.val     = {UGL };
noise.help    = {
                 'There are five noise prior options here (1) UGL, (2) GMRF, (3) LORETA '
                 '(4) Tissue-type and (5) Robust'
}';
noise.values  = {UGL GMRF LORETA tissue_type Robust };

%--------------------------------------------------------------------------
% first First level
%--------------------------------------------------------------------------
first         = cfg_menu;
first.tag     = 'first';
first.name    = 'First level';
first.val = {'No'};
first.help    = {
                 'This is implemented using Bayesian model comparison. For example, to test for the main effect of a factor two models are compared, one where the levels are represented using different regressors and one using the same regressor. This therefore requires explicit fitting of several models at each voxel and is computationally demanding (requiring several hours of computation). The recommended option is therefore NO.'
                 ''
                 'To use this option you must have already specified your factorial design during the model specification stage. '
}';
first.labels  = {'No', 'Yes'};
first.values  = {'No', 'Yes'};

%--------------------------------------------------------------------------
% second Second level
%--------------------------------------------------------------------------
second         = cfg_menu;
second.tag     = 'second';
second.name    = 'Second level';
second.val     = {'Yes'};
second.help    = {
                  'This option tells SPM to automatically generate the simple contrasts that are necessary to produce the contrast images for a second-level (between-subject) ANOVA. Naturally, these contrasts can also be used to characterise simple effects for each subject. '
                  ''
                  'With the Bayesian estimation option it is recommended that contrasts are computed during the parameter estimation stage (see ''simple contrasts'' below). The recommended option here is therefore YES.'
                  ''
                  'To use this option you must have already specified your factorial design during the model specification stage. '
                  ''
                  'If you wish to use these contrast images for a second-level analysis then you will need to spatially smooth them to take into account between-subject differences in functional anatomy ie. the fact that one persons V5 may be in a different position than anothers. '
}';
second.labels = {'No', 'Yes'};
second.values = {'No', 'Yes'};

%--------------------------------------------------------------------------
% anova ANOVA
%--------------------------------------------------------------------------
anova        = cfg_branch;
anova.tag    = 'anova';
anova.name   = 'ANOVA';
anova.val    = {first second };
anova.help   = {'Perform 1st or 2nd level Analysis of Variance.'};

%--------------------------------------------------------------------------
% name Name
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of contrast eg. ''Positive Effect'''};
name.strtype = 's';
name.num     = [1 Inf];

%--------------------------------------------------------------------------
% convec Contrast vector
%--------------------------------------------------------------------------
convec         = cfg_entry;
convec.tag     = 'convec';
convec.name    = 'Contrast vector';
convec.help    = {'These contrasts are used to generate PPMs which characterise effect sizes at each voxel. This is in contrast to SPMs in which eg. maps of t-statistics show the ratio of the effect size to effect variability (standard deviation). SPMs are therefore a-dimensional. This is not the case for PPMs as the size of the effect is of primary interest. Some care is therefore needed about the scaling of contrast vectors. For example, if you are interested in the differential effect size averaged over conditions then the contrast 0.5 0.5 -0.5 -0.5 would be more suitable than the 1 1 -1 -1 contrast which looks at the differential effect size summed over conditions. '};
convec.strtype = 'r';
convec.num     = [Inf 1];

%--------------------------------------------------------------------------
% gcon Simple contrast
%--------------------------------------------------------------------------
gcon         = cfg_branch;
gcon.tag     = 'gcon';
gcon.name    = 'Simple contrast';
gcon.val     = {name convec };
gcon.help    = {''};

%--------------------------------------------------------------------------
% generic Simple contrasts
%--------------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Simple contrasts';
generic.help    = {
                   '''Simple'' contrasts refers to a contrast that spans one-dimension ie. to assess an effect that is increasing or decreasing.'
                   ''
                   'If you have a factoral design then the contrasts needed to generate the contrast images for a 2nd-level ANOVA (or to assess these simple effects within-subject) can be specified automatically using the ANOVA->Second level option.'
                   ''
                   'When using the Bayesian estimation option it is computationally more efficient to compute the contrasts when the parameters are estimated. This is because estimated parameter vectors have potentially different posterior covariance matrices at different voxels and these matrices are not stored. If you compute contrasts post-hoc these matrices must be recomputed (an approximate reconstruction based on a Taylor series expansion is used). It is therefore recommended to specify as many contrasts as possible prior to parameter estimation.'
                   ''
                   'If you wish to use these contrast images for a second-level analysis then you will need to spatially smooth them to take into account between-subject differences in functional anatomy ie. the fact that one persons V5 may be in a different position than anothers. '
}';
generic.values  = {gcon };
generic.num     = [0 Inf];

%==========================================================================
% Bayesian Bayesian 1st-level
%==========================================================================
Bayesian         = cfg_branch;
Bayesian.tag     = 'Bayesian';
Bayesian.name    = 'Bayesian 1st-level';
Bayesian.val     = {space signal ARP noise LogEv anova generic};
Bayesian.help    = {
                    'Model parameters are estimated using Variational Bayes (VB). This allows you to specify spatial priors for regression coefficients and regularised voxel-wise AR(P) models for fMRI noise processes. The algorithm does not require functional images to be spatially smoothed. Estimation will take about 5 times longer than with the classical approach. This is why VB is not the default estimation option. '
                    ''
                    'Model estimation using this option is only efficient if MATLAB can load a whole slice of data into physical memory. With modern PCs this is usually achieved if the within-plane voxel sizes are 3 by 3 mm. This is therefore the minimum recommended voxel size that your spatial normalisation process should produce. Within-plane voxel sizes of 2 by 2 mm usually result in too many voxels per slice and result in estimation times lasting several hours or days. Such a high resolution is therefore to be avoided. '
                    ''
                    'After estimation, contrasts are used to find regions with effects larger than a user-specified size eg. 1 per cent of the global mean signal. These effects are assessed statistically using a Posterior Probability Map (PPM).'
}';

%==========================================================================
% Bayesian2 Bayesian 2nd-level
%==========================================================================
Bayesian2         = cfg_const;
Bayesian2.tag     = 'Bayesian2';
Bayesian2.name    = 'Bayesian 2nd-level';
Bayesian2.val     = {1};
Bayesian2.help    = {'Bayesian estimation of 2nd level models. This option uses the Empirical Bayes algorithm with global shrinkage priors that was previously implemented in SPM2. Use of the global shrinkage prior embodies a prior belief that, on average over all voxels, there is no net experimental effect. Some voxels will respond negatively and some positively with a variability determined by the prior precision. This prior precision can be estimated from the data using Empirical Bayes. '};

%==========================================================================
% method Method
%==========================================================================
method         = cfg_choice;
method.tag     = 'method';
method.name    = 'Method';
method.val     = {Classical};
method.help    = {
                  'There are three possible estimation procedures for fMRI models (1) classical (ReML) estimation of first or second level models, (2) Bayesian estimation of first level models and (3) Bayesian estimation of second level models. Option (2) uses spatial or global shrinkage priors. Option (3) uses global shrinkage priors. '
                  ''
                  'To use option (3) you must have already estimated the model using option (1). That is, for second-level models you must run a ReML estimation before running a Bayesian estimation. This is not necessary for option (2). Bayesian estimation of 1st-level models using VB does not require a prior ReML estimation.'
}';
method.values  = {Classical Bayesian Bayesian2};

%--------------------------------------------------------------------------
% write_residuals Write Residuals
%--------------------------------------------------------------------------
write_residuals        = cfg_menu;
write_residuals.tag    = 'write_residuals';
write_residuals.name   = 'Write residuals';
write_residuals.val    = {0};
write_residuals.help   = {'Write images of residuals to disk. This is only implemented for classical inference.'};
write_residuals.labels = {'No', 'Yes'};
write_residuals.values = {0, 1};

%--------------------------------------------------------------------------
% fmri_est Model estimation
%--------------------------------------------------------------------------
fmri_est          = cfg_exbranch;
fmri_est.tag      = 'fmri_est';
fmri_est.name     = 'Model estimation';
fmri_est.val      = {spmmat write_residuals method};
fmri_est.help     = {'Model parameters can be estimated using classical (ReML - Restricted Maximum Likelihood) or Bayesian algorithms. After parameter estimation, the RESULTS button can be used to specify contrasts that will produce Statistical Parametric Maps (SPMs) or Posterior Probability Maps (PPMs) and tables of statistics.'};
fmri_est.prog     = @spm_run_fmri_est;
fmri_est.vout     = @vout_stats;
fmri_est.modality = {'FMRI' 'PET' 'EEG'};


%==========================================================================
% function dep = vout_stats(job)
%==========================================================================
function dep = vout_stats(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'SPM.mat File';
dep(1).src_output = substruct('.','spmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%dep(2)            = cfg_dep;
%dep(2).sname      = 'SPM Variable';
%dep(2).src_output = substruct('.','spmvar');
%dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
if isfield(job.method, 'Classical')
    dep(2)            = cfg_dep;
    dep(2).sname      = 'Beta Images';
    dep(2).src_output = substruct('.','beta');
    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(3)            = cfg_dep;
    dep(3).sname      = 'Analysis Mask';
    dep(3).src_output = substruct('.','mask');
    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(4)            = cfg_dep;
    dep(4).sname      = 'ResMS Image';
    dep(4).src_output = substruct('.','resms');
    dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if job.write_residuals
        dep(4)            = cfg_dep;
        dep(4).sname      = 'Residual Images';
        dep(4).src_output = substruct('.','res');
        dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    % can't check whether auto-generated contrasts are generated this is
    % specified in input SPM.mat, not this job
end
