function fmri = spm_cfg_dcm_fmri
% SPM Configuration file for DCM for fMRI
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: spm_cfg_dcm_fmri.m 6711 2016-02-03 15:25:43Z peter $

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

%--------------------------------------------------------------------------
% dir Directory
%--------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select the directory where the output will be written.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

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
% subj Create single model
%--------------------------------------------------------------------------
model      = cfg_branch;
model.tag  = 'model';
model.name = 'Model';
model.val  = {dcmmat};
model.help = {'Corresponding model for each subject'};

% -------------------------------------------------------------------------
% subjects Create set of models
%--------------------------------------------------------------------------
models        = cfg_repeat;
models.tag    = 'models';
models.name   = 'Per model';
models.values = {model};
models.help   = {'Select DCM.mat files per model'};
models.num    = [1 Inf];

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
% fmri Dynamic Causal Model for fMRI
% -------------------------------------------------------------------------
fmri         = cfg_choice; 
fmri.tag     = 'fmri';
fmri.name    = 'DCM for fMRI';
fmri.help    = {'Dynamic Causal Modelling for fMRI'};
fmri.values  = { regions inputs };


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