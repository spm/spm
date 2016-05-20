function cfg = tbx_cfg_longitudinal
% MATLABBATCH Configuration file for toolbox 'Longitudinal'
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_longitudinal.m 6798 2016-05-20 11:53:33Z john $

if ~isdeployed,
    addpath(fullfile(spm('Dir'),'toolbox','Longitudinal'));
    addpath(fullfile(spm('dir'),'toolbox','Shoot'));
end

% ---------------------------------------------------------------------
% vols Volumes
% ---------------------------------------------------------------------
vols1         = cfg_files;
vols1.tag     = 'vols1';
vols1.name    = 'Time 1 Volumes';
vols1.help    = {'Select first time point scans of each subject.'};
vols1.filter = 'image';
vols1.ufilter = '.*';
vols1.num     = [1 Inf];

vols2         = cfg_files;
vols2.tag     = 'vols2';
vols2.name    = 'Time 2 Volumes';
vols2.help    = {'Select second time point scans of each subject. Note that the order that the first and second time points are specified should be the same.  The algorithm does not incorporate any magical way of figuring out which scans go together.'};
vols2.filter = 'image';
vols2.ufilter = '.*';
vols2.num     = [1 Inf];

vols         = cfg_files;
vols.tag     = 'vols';
vols.name    = 'Volumes';
vols.help    = {'Select scans for this subject.'};
vols.filter = 'image';
vols.ufilter = '.*';
vols.num     = [1 Inf];


tdif         = cfg_entry;
tdif.tag     = 'tdif';
tdif.name    = 'Time Difference';
tdif.help    = {'Specify the time difference between the scans in years.  This can be a single value (if it is the same for all subjects) or a vector of values (if it differs among subjects).'};
tdif.strtype = 'e';
tdif.num     = [1 Inf];
tdif.val     = {1};

tim         = cfg_entry;
tim.tag     = 'times';
tim.name    = 'Times';
tim.help    = {'Specify the times of the scans in years.'};
tim.strtype = 'e';
tim.num     = [1 Inf];

noise         = cfg_entry;
noise.tag     = 'noise';
noise.name    = 'Noise Estimate';
noise.help    = {'.'};
noise.strtype = 'e';
noise.num     = [Inf Inf];
noise.val     = {NaN};
noise.help    = {'Specify the standard deviation of the noise in the images.  If a scalar is entered, all images will be assumed to have the same level of noise.  For any non-finite values, the algorithm will try to estimate reasonable noise estimates based on fitting a mixture of two Rician distributions to the intensity histogram of each of the images. This works reasonably well for simple MRI scans, but less well for derived images (such as averages).  The assumption is that the residuals, after fitting the registration model, are i.i.d. Gaussian.'
};

bparam         = cfg_entry;
bparam.tag     = 'bparam';
bparam.name    = 'Bias Regularisation';
bparam.help    = {'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   ''
                   'An important issue relates to the distinction between variations in the difference between the images that arise because of the differential bias artifact due to the physics of MR scanning, and those that arise due to shape differences.  The objective is to model the latter by deformations, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large estimates of the intensity non-uniformity.'
                   'Knowing what works best should be a matter of empirical exploration, as it depends on the scans themselves.  For example, if your data has very little of the artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
                   }';
bparam.strtype = 'e';
bparam.num      = [1 1];
bparam.val      = {1e6};

wparam         = cfg_entry;   
wparam.tag     = 'wparam';
wparam.name    = 'Warping Regularisation';
wparam.help    = {'Registration involves simultaneously minimising two terms.  One of these is a measure of similarity between the images (mean-squared difference in the current situation), whereas the other is a measure of the roughness of the deformations.  This measure of roughness involves the sum of the following terms:',...
'* Absolute displacements need to be penalised by a tiny amount.  The first element encodes the amount of penalty on these.  Ideally, absolute displacements should not be penalised, but it is often necessary for technical reasons.',...
'* The `membrane energy'' of the deformation is penalised (2nd element), usually by a relatively small amount. This penalises the sum of squares of the derivatives of the velocity field (ie the sum of squares of the elements of the Jacobian tensors).',...
'* The `bending energy'' is penalised (3rd element). This penalises the sum of squares of the 2nd derivatives of the velocity.',...
'* Linear elasticity regularisation is also included (4th and 5th elements).  The first parameter (mu) is similar to that for linear elasticity, except it penalises the sum of squares of the Jacobian tensors after they have been made symmetric (by averaging with the transpose).  This term essentially penalises length changes, without penalising rotations.',...
'* The final term also relates to linear elasticity, and is the weight that denotes how much to penalise changes to the divergence of the velocities (lambda).  This divergence is a measure of the rate of volumetric expansion or contraction.',...
'Note that regularisation is specified based on what is believed to be appropriate for a year of growth.  The specified values are divided by the number of years time difference.' 
};
wparam.strtype = 'e';
wparam.num     = [1 5];
wparam.val     = {[0 0 100 25 100]};
% Change to (eg): wparam.val     = {[0 0 100 25 12]};

write_avg         = cfg_menu;
write_avg.tag     = 'write_avg';
write_avg.name    = 'Save Mid-point average';
write_avg.help    = {'Do you want to save the mid-point average template image? This is likely to be useful for groupwise alignment, and is prefixed by ``avg_'''' and written out in the same directory of the first time point data.'};
write_avg.labels = {
                'Save'
                'Dont save'
                }';
write_avg.values = { 1 0 };
write_avg.val    = {1};

write_div         = cfg_menu;
write_div.tag     = 'write_div';
write_div.name    = 'Save Divergence';
write_div.help    = {'Do you want to save a map of divergence of the velocity field?  This is useful for morphometrics, and may be considered as the rate of volumetric expansion.  Negative values indicate contraction. These files are prefixed by ``dv_'''' and written out in the same directory of the first time point data.'};
write_div.labels = {
                'Save'
                'Dont save'
                }';
write_div.values = { 1 0 };
write_div.val    = {1};

write_divd         = write_div;
write_divd.name    = 'Save Divergence Rate';
write_divd.help    = {'Do you want to save a map of divergence of the velocity field?  This is useful for morphometrics, and may be considered as the rate of volumetric expansion.  Negative values indicate contraction. These files are prefixed by ``dv_'''' and written out in the same directory of the first time point data. Note that the divergences written out have been divided by the time interval between scans'};

write_jacd         = cfg_menu;
write_jacd.tag     = 'write_jac';
write_jacd.name    = 'Save Jacobian Rate';
write_jacd.help    = {'Do you want to save a map of the differences between the Jacobian determinants, divided by the time interval?  Some consider these useful for morphometrics (although the divergences of the initial velocities may be preferable). The difference between two Jacobian determinants is computed and this is divided by the time interval. One original Jacobian map is for the deformation from the mid point to the first scan, and the other is for the deformation from the mid point to the second scan.  Each of these encodes the relative volume (at each spatial location) between the scan and the mid-point average. Values less than 0 indicate contraction (over time), whereas values greater than zero indicate expansion.  These files are prefixed by ``jd_'''' and written out in the same directory of the first time point data.'};
write_jacd.labels = {
                'Save'
                'Dont save'
                }';
write_jacd.values = { 1 0 };
write_jacd.val    = {0};

write_jac         = cfg_menu;
write_jac.tag     = 'write_jac';
write_jac.name    = 'Save Jacobians';
write_jac.help    = {'Do you want to save a map of the Jacobian determinants?  Some consider these useful for morphometrics (although the divergences of the initial velocities may be preferable). Each map of Jacobians encodes the relative volume (at each spatial location) between the scan and the median time-point average. Values less than one indicate contraction (over time), whereas values greater than one indicate expansion.  These files are prefixed by ``j_'''' and written out in the same directory of the first time point data.'};
write_jac.labels = {
                'Save'
                'Dont save'
                }';
write_jac.values = { 1 0 };
write_jac.val    = {0};

write_defs         = cfg_menu;
write_defs.tag     = 'write_def';
write_defs.name    = 'Deformation Fields';
write_defs.help    = {'Deformation fields can be saved to disk, and used by the Deformations Utility. Deformations are saved as y_*.nii files, which contain three volumes to encode the x, y and z coordinates.  They are written in the same directory as the corresponding image.'};
write_defs.labels = {
                'Save'
                'Dont save'
                }';
write_defs.values = { 1 0 };
write_defs.val    = {0};

% ---------------------------------------------------------------------
% longit2 Pairwise Longitudinal Registration
% ---------------------------------------------------------------------
long2         = cfg_exbranch;
long2.tag     = 'pairwise';
long2.name    = 'Pairwise Longitudinal Registration';
long2.val     = {vols1 vols2 tdif noise wparam bparam write_avg write_jacd write_divd write_defs};
long2.help    = {'Longitudinal registration of pairs of anatomical MRI scans.  It is based on pairwise inverse-consistent alignment between the first and second scan of each subject, and incorporates a bias field correction.  Prior to running the registration, the scans should already be in very rough alignment, although because the model incorporates a rigid-body transform, this need not be extremely precise.  Note that there are a bunch of hyper-parameters to be specified.  If you are unsure what values to take, then the defaults should be a reasonable guess of what works.  Note that changes to these hyper-parameters will impact the results obtained.',...
'',...
'The alignment assumes that all scans have similar resolutions and dimensions, and were collected on the same (or very similar) MR scanner using the same pulse sequence.  If these assumption are not correct, then the approach will not work as well.'};
long2.prog = @spm_pairwise;
long2.vout = @vout_pairwise;
%----------------------------------------------------------------------

long         = cfg_exbranch;
long.tag     = 'series';
long.name    = 'Serial Longitudinal Registration';
long.val     = {vols tim noise wparam bparam write_avg write_jac write_div write_defs};
long.help    = {'Longitudinal registration of series of anatomical MRI scans for a single subject.  It is based on groupwise alignment among each of the subject''s scans, and incorporates a bias field correction.  Prior to running the registration, the scans should already be in very rough alignment, although because the model incorporates a rigid-body transform, this need not be extremely precise.  Note that there are a bunch of hyper-parameters to be specified.  If you are unsure what values to take, then the defaults should be a reasonable guess of what works.  Note that changes to these hyper-parameters will impact the results obtained.',...
'',...
'The alignment assumes that all scans have similar resolutions and dimensions, and were collected on the same (or very similar) MR scanner using the same pulse sequence.  If these assumption are not correct, then the approach will not work as well.  There are a number of settings (noise estimate, regularisation etc). Default settings often work well, but it can be very helpful to try some different values, as these can have a large effect on the results.'};
long.prog = @spm_series_align;

cfg        = cfg_repeat;
cfg.tag    = 'longit';
cfg.name   = 'Longitudinal Registration';
cfg.values = {long2,long};


%======================================================================
function cdep = vout_pairwise(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

ind  = 1;
if job.write_avg,
    cdep(ind)          = cfg_dep;
    cdep(ind).sname      = 'Midpoint Average';
    cdep(ind).src_output = substruct('.','avg','()',{':'});
    cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    ind = ind + 1;
end
if job.write_jac,
    cdep(ind)          = cfg_dep;
    cdep(ind).sname      = 'Jacobian Diff';
    cdep(ind).src_output = substruct('.','jac','()',{':'});
    cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    ind = ind + 1;
end
if job.write_div,
    cdep(ind)          = cfg_dep;
    cdep(ind).sname      = 'Divergence';
    cdep(ind).src_output = substruct('.','div','()',{':'});
    cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    ind = ind + 1;
end
if job.write_def,
    cdep(ind)            = cfg_dep;
    cdep(ind).sname      = 'Deformation (1)';
    cdep(ind).src_output = substruct('.','def1','()',{':'});
    cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    ind = ind + 1;

    cdep(ind)            = cfg_dep;
    cdep(ind).sname      = 'Deformation (2)';
    cdep(ind).src_output = substruct('.','def2','()',{':'});
    cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    ind = ind + 1;
end

