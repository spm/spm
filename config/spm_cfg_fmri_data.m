function fmri_data = spm_cfg_fmri_data
% SPM Configuration file for fMRI data specification
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_fmri_data.m 6088 2014-07-03 17:57:09Z guillaume $


%--------------------------------------------------------------------------
% scans Scans
%--------------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the fMRI scans for this session.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter  = {'image','mesh'};
scans.ufilter = '.*';
scans.num     = [1 Inf];

%--------------------------------------------------------------------------
% spmmat Select SPM.mat
%--------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file containing the specified design matrix.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '.*';
spmmat.num     = [1 1];

%--------------------------------------------------------------------------
% mask Explicit mask
%--------------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Explicit mask';
mask.val{1}  = {''};
mask.help    = {'Specify an image for explicitly masking the analysis. A sensible option here is to use a segmention of structural images to specify a within-brain mask. If you select that image as an explicit mask then only those voxels in the brain will be analysed. This both speeds the estimation and restricts SPMs/PPMs to within-brain voxels. Alternatively, if such structural images are unavailable or no masking is required, then leave this field empty.'};
mask.filter  = {'image','mesh'};
mask.ufilter = '.*';
mask.num     = [0 1];

%--------------------------------------------------------------------------
% fmri_data fMRI data specification
%--------------------------------------------------------------------------
fmri_data          = cfg_exbranch;
fmri_data.tag      = 'fmri_data';
fmri_data.name     = 'fMRI data specification';
fmri_data.val      = {scans spmmat mask };
fmri_data.help     = {'Select the data and optional explicit mask for a specified design'};
fmri_data.prog     = @spm_run_fmri_data;
fmri_data.vout     = @vout_stats;
fmri_data.modality = {'FMRI'};


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
