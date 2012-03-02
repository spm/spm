function long = tbx_cfg_longitudinal
% MATLABBATCH Configuration file for toolbox 'Longitudinal'
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_longitudinal.m 4671 2012-03-02 19:40:35Z john $

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','Longitudinal')); end

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
vols2.help    = {'Select second time point scans of each subject.'};
vols2.filter = 'image';
vols2.ufilter = '.*';
vols2.num     = [1 Inf];

tdif         = cfg_entry;
tdif.tag     = 'tdif';
tdif.name    = 'Time difference';
tdif.help    = {'Specify the time difference between the scans in years.  This can be a single value (if it is the same for all subjects) or a vector of values (if it differs among subjects).'};
tdif.strtype = 'e';
tdif.num     = [1 Inf];
tdif.val     = {1};

noise         = cfg_entry;
noise.tag     = 'noise';
noise.name    = 'Noise Estimate';
noise.help    = {'.'};
noise.strtype = 'e';
noise.num     = [1 1];
noise.val     = {NaN};

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
wparam.help    = {'Registration involves simultaneously minimising two terms.  One of these is a measure of similarity between the images (mean-squared difference in the current situation), whereas the other is a measure of the roughness of the deformations.  This measure of roughness involves the sum of the following terms:'
'* Absolute displacements need to be penalised by a tiny amount.  The first element encodes the amount of penalty on these.  Ideally, absolute displacements should not be penalised, but it is necessary for technical reasons.'
'* The `membrane energy'' of the deformation is penalised (2nd element), usually by a relatively small amount. This penalises the sum of squares of the derivatives of the velocity field (ie the sum of squares of the elements of the Jacobian tensors).'
'* The `bending energy'' is penalised (3rd element). This penalises the sum of squares of the 2nd derivatives of the velocity.'
'* Linear elasticity regularisation is also included (4th and 5th elements).  The first parameter (mu) is similar to that for linear elasticity, except it penalises the sum of squares of the Jacobian tensors after they have been made symmetric (by averaging with the transpose).  This term essentially penalises length changes, without penalising rotations.'
'* The final term also relates to linear elasticity, and is the weight that denotes how much to penalise changes to the divergence of the velocities (lambda).  This divergence is a measure of the rate of volumetric expansion or contraction.' 
};
wparam.strtype = 'e';
wparam.num     = [1 5];
wparam.val     = {[0.05 5 100 50 200]};

write_avg         = cfg_menu;
write_avg.tag     = 'write_avg';
write_avg.name    = 'Save Mid-point average';
write_avg.help    = {'Do you want to save the mid-point average template image? This is likely to be useful for groupwise alignment.'};
write_avg.labels = {
                'Save average'
                'Dont save average'
                }';
write_avg.values = {
                [1]
                [0]
                }';
write_avg.val    = {[1]};

write_div         = cfg_menu;
write_div.tag     = 'write_div';
write_div.name    = 'Save divergence';
write_div.help    = {'Do you want to save a map of divergence of the velocity field?  This is useful for morphometrics, and may be considered as the rate of volumetric expansion.  Negative values indicate contraction.'};
write_div.labels = {
                'Save divergence'
                'Dont save divergence'
                }';
write_div.values = {
                [1]
                [0]
                }';
write_div.val    = {[1]};

write_jac         = cfg_menu;
write_jac.tag     = 'write_jac';
write_jac.name    = 'Save Jacobians';
write_jac.help    = {'Do you want to save a map of Jacobian determinants?  This is useful for morphometrics, and may be considered as the relative volume at each point in the two images.  Values less than 1 indicate contraction.'};
write_jac.labels = {
                'Save Jacobians'
                'Dont save Jacobians'
                }';
write_jac.values = {
                [1]
                [0]
                }';
write_jac.val    = {[1]};

write_defs         = cfg_menu;
write_defs.tag     = 'write_def';
write_defs.name    = 'Deformation Fields';
write_defs.help    = {'Deformation fields can be saved to disk, and used by the Deformations Utility. Deformations are saved as y_*.nii files, which contain three volumes to encode the x, y and z coordinates.'};
write_defs.labels = {
                'Dont save deformations'
                'Save deformations'
                }';
write_defs.values = {
                [0]
                [1]
                }';
write_defs.val    = {[0]};

%write         = cfg_menu;
%write.tag     = 'write';
%write.name    = 'Deformation Fields';
%write.help    = {'Deformation fields can be saved to disk, and used by the Deformations Utility. For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'};
%write.labels = {
%                'None'
%                'Time 2 -> Time 1'
%                'Time 1 -> Time 2'
%                'Both'
%                }';
%write.values = {
%                [0 0]
%                [1 0]
%                [0 1]
%                [1 1]
%                }';
%write.val    = {[0 0]};
% ---------------------------------------------------------------------
% long Longitudinal Registration
% ---------------------------------------------------------------------
long         = cfg_exbranch;
long.tag     = 'longit2';
long.name    = 'Longitudinal Registration';
long.val     = {vols1 vols2 tdif noise wparam bparam write_avg write_jac write_div write_defs};
long.help    = {'This toolbox is currently work in progress.'};
long.prog = @spm_pairwise;
long.vout = @vout;
%----------------------------------------------------------------------


%======================================================================
function cdep = vout(job)
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
    cdep(ind).sname      = 'Jacobian';
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

