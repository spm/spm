function setlevel = spm_cfg_setlevel
% SPM Configuration file for Set level tests based on Barnes et al. NIMG
% 2012
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_setlevel.m 5551 2013-06-12 13:14:04Z gareth $

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

% ---------------------------------------------------------------------
% Set level test 
% ---------------------------------------------------------------------
setlevel          = cfg_exbranch;
setlevel.tag      = 'setlevel';
setlevel.name     = 'Set Level test';
setlevel.val      = {spmmat};
setlevel.help     = {'A set level test how likely the statistical image is a random field'};
setlevel.prog     = @spm_run_setlevel;
%setlevel.vout     = @vout_stats;
setlevel.modality = {'FMRI' 'PET' 'EEG'};


