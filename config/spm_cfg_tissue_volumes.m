function tvol = spm_cfg_tissue_volumes
% SPM Configuration file for Tissue Volumes
%
% See also: spm_run_tissue_volumes, spm_summarise
%__________________________________________________________________________
% Copyright (C) 2013-2016 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_cfg_tissue_volumes.m 6952 2016-11-25 16:03:13Z guillaume $


mat         = cfg_files;
mat.tag     = 'matfiles';
mat.name    = 'Segmentation mat-files';
mat.filter  = 'mat';
mat.ufilter = 'seg8\.mat$';
mat.num     = [1 Inf];
mat.help    = {
    'Select the ''seg8.mat'' files containing the segmentation parameters.'
    };

T           = cfg_entry;
T.tag       = 'tmax';
T.name      = 'Maximum tissue class';
T.strtype   = 'n';
T.num       = [1 1];
T.val       = {3};
T.help      = {
    ['Specify the maximum tissue class, T, where tissues 1:T will be ' ...
    'measured.']
    ['The default of 3 corresponds to GM, WM and CSF for the ' ...
    'default tissue prior probability maps ''TPM.nii,1'' to ''TPM.nii,3''']
    ''
    ['The sum of these tissues will also be computed, which by ' ...
    'default is the total intracranial volume (known as TIV or ICV). ' ...
    'If T=2, the sum will by default be the total parenchymal brain ' ...
    'volume (known as TBV or PBV), which is also often of interest.']
    };

mask        = cfg_files;
mask.tag    = 'mask';
mask.name   = 'Mask image';
mask.filter = 'image';
mask.num    = [0 1];
mask.val    = {{fullfile(spm('dir'), 'tpm', 'mask_ICV.nii,1')}};
mask.help   = {
    'Optional binary mask image in same space as TPMs (e.g. MNI).'
    'Only voxels inside this mask will count for the total volume.'
    ''
    ['The default mask excludes the eyes, which might otherwise be ' ...
    'counted in the fluid tissue class (that includes cerebrospinal, ' ...
    'aqueous, and vitreous fluids/humours).']
    };

outf         = cfg_entry;
outf.tag     = 'outf';
outf.name    = 'Output file';
outf.strtype = 's';
outf.val     = {''};
outf.help    = {
    'Filename for saving results.'
    ['Segmentation filenames and volumes ' ...
    'will be stored in CSV format (comma-separated variables).']
    ''
    'This can be empty; results will appear in the MATLAB command window.'
    };

tvol        = cfg_exbranch;
tvol.tag    = 'tvol';
tvol.name   = 'Tissue Volumes';
tvol.val    = {mat T mask outf};
tvol.help   = {
    'Compute total tissue volumes (in litres) from segmentation results.'
    ''
    ['Only the seg8.mat files are required, but if modulated warped ' ...
    'segmentations (mwc*) are found they will be reused, saving time ' ...
    '(and allowing you to use non-default values for MRF and/or ' ...
    'clean-up options if you wish).']
    };
tvol.prog = @(job) spm_run_tissue_volumes('exec', job);
tvol.vout = @(job) spm_run_tissue_volumes('vout', job);
