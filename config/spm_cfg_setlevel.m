function setlevel = spm_cfg_setlevel
% SPM Configuration file for Set level tests based on Barnes et al. NIMG
% 2012
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_setlevel.m 5554 2013-06-13 09:13:56Z gareth $

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {
                  'Select the SPM.mat file that contains the design matrix specification and results. '
                  
}';
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];


cindex = cfg_entry;
cindex.tag = 'cindex';
cindex.name = 'Contrast index';
cindex.strtype = 'n';
cindex.help = {'Index of the contrast of interest'};
cindex.val = {1};


% ---------------------------------------------------------------------
% Set level test 
% ---------------------------------------------------------------------
setlevel          = cfg_exbranch;
setlevel.tag      = 'setlevel';
setlevel.name     = 'Set Level test';
setlevel.val      = {spmmat cindex};
setlevel.help     = {'A set level test how likely the statistical image is a random field'};
setlevel.prog     = @spm_run_setlevel;
%setlevel.vout     = @vout_stats;
setlevel.modality = {'FMRI' 'PET' 'EEG'};


