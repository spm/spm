function fmri = spm_cfg_dcm_fmri
% SPM Configuration file for DCM for fMRI
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_dcm_fmri.m 6565 2015-09-30 10:42:14Z peter $

% -------------------------------------------------------------------------
% dcmmat Select DCM_*.mat
% -------------------------------------------------------------------------
dcmmat         = cfg_files;
dcmmat.tag     = 'dcmmat';
dcmmat.name    = 'Select DCM_*.mat';
dcmmat.help    = {'Select DCM_*.mat files.'};
dcmmat.filter  = 'mat';
dcmmat.ufilter = '^DCM_.*\.mat$';
dcmmat.num     = [1 Inf];

% -------------------------------------------------------------------------
% voimat Select VOI_*.mat
% -------------------------------------------------------------------------
voimat         = cfg_files;
voimat.tag     = 'voimat';
voimat.name    = 'Select VOI_*.mat';
voimat.help    = {'Select VOI_*.mat files.'};
voimat.filter  = 'mat';
voimat.ufilter = '^VOI_.*\.mat$';
voimat.num     = [1 Inf];

% -------------------------------------------------------------------------
% spmmat Select SPM.mat
% -------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select SPM.mat file.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% -------------------------------------------------------------------------
% session Session index
% -------------------------------------------------------------------------
session         = cfg_entry;
session.tag     = 'session';
session.name    = 'Which session';
session.help    = {'Enter the session number.'};
session.strtype = 'e';
session.num     = [1 1];

% -------------------------------------------------------------------------
% val Val
% -------------------------------------------------------------------------
val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Values';
val.help    = {'Inputs to include for one condition. Enter ''1'' ' ...
               'to include this condition (with no parameteric regressor). '...
               'Entering [1 0 1] would include this condition and '...
               'its second parametric regressor.'};
val.strtype = 'e';
val.num     = [1 Inf];

% -------------------------------------------------------------------------
% inp Inputs
% -------------------------------------------------------------------------
inp         = cfg_repeat;
inp.tag     = 'inputs';
inp.name    = 'Inputs';
inp.help    = {'Inputs to include and their parametric modulations (PMs). '...
               'You should click ''New: Values'' for each condition in '...
               'your SPM (i.e. SPM.U), up to the last condition you wish  '...
               'to include.'};
inp.values  = { val };
inp.num     = [1 Inf];

% -------------------------------------------------------------------------
% regions Specify regions
% -------------------------------------------------------------------------
regions      = cfg_exbranch;
regions.tag  = 'regions';
regions.name = 'Region specification';
regions.val  = { dcmmat voimat };
regions.help = {'Insert new regions into a DCM model. '...
    '' ...
    'The RT is assumed to be the same as before. '...
    ''...
    ['This functionality can be used, for example, to replace subject X''s '...
    'data by subject Y''s. The model can then be re-estimated without '...
    'having to go through model specification again.']};
regions.prog = @spm_run_dcm_fmri_regions;
regions.vout = @vout_dcm_fmri;

% -------------------------------------------------------------------------
% inputs Specify inputs
% -------------------------------------------------------------------------
inputs      = cfg_exbranch;
inputs.tag  = 'inputs';
inputs.name = 'Input specification';
inputs.val  = { dcmmat spmmat session inp };
inputs.help = {'Insert new inputs into a DCM model'...
    ''...
    ['This functionality can be used, for example, to replace subject X''s '...
    'inputs by subject Y''s. The model can then be re-estimated without '...
    'having to go through model specification again.']};
inputs.prog = @spm_run_dcm_fmri_inputs;
inputs.vout = @vout_dcm_fmri;

% -------------------------------------------------------------------------
% analysis Analysis
% -------------------------------------------------------------------------
analysis         = cfg_menu;
analysis.tag     = 'analysis';
analysis.name    = 'Analysis';
analysis.help    = {'Analysis model.'};
analysis.labels  = {'time series','cross-spectral densities'};
analysis.values  = {'time','csd'};
analysis.val     = {'time'};

% -------------------------------------------------------------------------
% estimate Estimate
% -------------------------------------------------------------------------
estimate      = cfg_exbranch;
estimate.tag  = 'estimate';
estimate.name = 'DCM estimation';
estimate.val  = { dcmmat analysis };
estimate.help = {'Estimate parameters of a DCM for fMRI data.'};
estimate.prog = @spm_run_dcm_fmri_est;
estimate.vout = @vout_dcm_fmri;

% -------------------------------------------------------------------------
% fmri Dynamic Causal Model for fMRI
% -------------------------------------------------------------------------
fmri         = cfg_choice; 
fmri.tag     = 'fmri';
fmri.name    = 'DCM for fMRI';
fmri.help    = {'Dynamic Causal Modelling for fMRI'};
fmri.values  = { regions inputs estimate };

%==========================================================================
function out = spm_run_dcm_fmri_est(job)
%==========================================================================
for i=1:numel(job.dcmmat)
    switch lower(job.analysis)
        case 'time'
            spm_dcm_estimate(job.dcmmat{i});
        case 'csd'
            spm_dcm_fmri_csd(job.dcmmat{i});
    end
end
out = job.dcmmat;

%==========================================================================
function out = spm_run_dcm_fmri_inputs(job)
%==========================================================================
for i=1:numel(job.dcmmat)
    spm_dcm_U(job.dcmmat{i},job.spmmat{1},job.session,job.val);
end
out = job.dcmmat;

%==========================================================================
function out = spm_run_dcm_fmri_regions(job)
%==========================================================================
for i=1:numel(job.dcmmat)
    spm_dcm_voi(job.dcmmat{i},job.voimat);
end
out = job.dcmmat;

%==========================================================================
function dep = vout_dcm_fmri(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'DCM mat File(s)';
dep(1).src_output = substruct('.','dcmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
