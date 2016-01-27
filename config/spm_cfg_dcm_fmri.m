function fmri = spm_cfg_dcm_fmri
% SPM Configuration file for DCM for fMRI
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: spm_cfg_dcm_fmri.m 6695 2016-01-27 10:51:26Z peter $

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
% dcmmat Select GCM_.*.mat
% -------------------------------------------------------------------------
gcmmat         = cfg_files;
gcmmat.tag     = 'gcmmat';
gcmmat.name    = 'Select GCM_*.mat';
gcmmat.help    = {'Select GCM_*.mat files.'};
gcmmat.filter  = 'mat';
gcmmat.ufilter = '^GCM_.*\.mat$';
gcmmat.num     = [1 1];

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
% name Model name
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {['Specify a name for the group DCM file. The prefix GCM_' ...
                'and suffix .mat are automatically added']};
name.strtype = 's';
name.num     = [0 Inf];

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
% subj Create single subject
%--------------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {dcmmat};
subj.help = {'Subject with one or more models.'};

% -------------------------------------------------------------------------
% multiple_models Create set of subjects
%--------------------------------------------------------------------------
subjects        = cfg_repeat;
subjects.tag    = 'subjects';
subjects.name   = 'Select multiple models per subject';
subjects.values = {subj};
subjects.help   = {'Create the subjects and select the models for each'};
subjects.num    = [1 Inf];

% -------------------------------------------------------------------------
% models Select models (one per subject)
%--------------------------------------------------------------------------
models      = cfg_branch;
models.tag  = 'models';
models.name = 'Select one model per subject';
models.val  = {dcmmat};
models.help = {'Models - one per subject.'};

% -------------------------------------------------------------------------
% output_single Output one .mat file for the group
%--------------------------------------------------------------------------
output_single         = cfg_branch;
output_single.tag     = 'output_single';
output_single.name    = 'Create group GCM_*.mat file';
output_single.val     = { dir name };
output_single.help    = {['Creates a single group-level DCM file ' ...
                          'containing a subjects x models cell array.']};

% -------------------------------------------------------------------------
% output_separate Output one .mat file per model
%--------------------------------------------------------------------------
output_separate         = cfg_branch;
output_separate.tag     = 'output_separate';
output_separate.name    = 'Output individual DCM files';
output_separate.val     = {};
output_separate.help    = {'Updated existing individual DCM.mat files'};

% -------------------------------------------------------------------------
% output_type Choice of how many DCM.mat files to output
%--------------------------------------------------------------------------
output_type         = cfg_choice;
output_type.tag     = 'output_type';
output_type.name    = 'Output';
output_type.values  = { output_single output_separate };
output_type.val     = { output_single };
output_type.help    = {['Whether to create a single DCM file across all ' ...
                       'subjects / models (default, required for second ' ...
                       'level analysis) or just update the separate' ...
                       'first-level DCM files.']};

% -------------------------------------------------------------------------
% way Choice of ways to select DCMs (nested models)
%--------------------------------------------------------------------------
dcms        = cfg_choice;
dcms.tag    = 'dcms';
dcms.name   = 'Select DCMs';
dcms.values = {models subjects gcmmat};
dcms.val    = {subjects};
dcms.help   = {['Select one DCM per subject, multiple DCMs per subject ' ...
                  'or an existing group DCM file.'] ...
                 ['If multiple DCMs are selected per subject, then ' ...
                  'the first DCM for each subject should a ''full'' ' ...
                  'model containing all connections of interest. Subsequent ' ...
                  '(nested) DCMs will have certain connections switched ' ...
                  'off.']};                   

% -------------------------------------------------------------------------
% est_type Estimation type
%--------------------------------------------------------------------------     
est_type        = cfg_menu;
est_type.tag    = 'est_type';
est_type.name   = 'Estimation type';
est_type.labels = {'Full + BMR (default)',...
                   'Full + BMR PEB (more accurate but slower)',...
                   'Full (not recommended)',...
                   'None (collate only)'};
est_type.values = {1,2,3,4};
est_type.val    = {1};
est_type.help   = {['Full + BMR: Estimates the full (first) model for ' ...
                    'each subject then uses Bayesian Model Reduction (BMR) '...
                    'to rapidly infer the evidence / parameters for any ' ...
                    'subsequent nested models.'] ...
                   ['Full + BMR PEB: Iteratively estimates the full (first) '...
                    'model for each subject, then sets the priors on each' ...
                    'each parameter to the group mean (from a PEB model)' ...
                    'then re-estimates. This improves estimation '...
                    'by overcoming local optima, but takes longer.' ] ...
                   ['Full: Estimates all models individually. Provided for '...
                    'backward compatibility.']...
                   ['None: Creates a group level DCM file without ' ...
                    'performing estimation']};
                              
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
analysis.labels  = {'time series','cross-spectral densities'};
analysis.values  = {'time','csd'};
analysis.val     = {'time'};
analysis.help    = {['Whether to analyse in the time domain (for task-' ...
                     'based studies or stochastic DCM) in the frequency '...
                     'domain (for resting state analysis with DCM for CSD']};

% -------------------------------------------------------------------------
% estimate Estimate
% -------------------------------------------------------------------------
estimate      = cfg_exbranch;
estimate.tag  = 'estimate';
estimate.name = 'DCM estimation';
estimate.val  = { output_type dcms analysis est_type };
estimate.help = {['Estimate the parameters and free energy (log model ' ...
                  'evidence) of first level DCMs for fMRI. Models ' ...
                  'are assembled into a Subjects x Models array and ' ...
                  'saved in group GCM_*.mat file']};
estimate.prog = @spm_run_dcm_fmri_est;
estimate.vout = @vout_gcm_fmri;

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

EST_FULL_BMR     = 1;
EST_FULL_BMR_PEB = 2;
EST_FULL         = 3;
EST_NONE         = 4;

dcms = job.dcms;

% Build subjects x models filename matrix
if isfield(dcms,'models')
    P = dcms.models.dcmmat;
    ns = size(P,1);
    nm = 1;
        
    % Load all models into memory
    GCM = spm_dcm_load(P);    
    
elseif isfield(dcms,'subj')
    ns  = length(dcms.subj);
    nm  = length(dcms.subj(1).dcmmat);
    P = cell(ns,nm);
    
    for s = 1:ns
        if length(dcms.subj(s).dcmmat) ~= nm
            error(['Please ensure all subjects have the same number of ' ... 
                   'models']);
        end
        
        P(s,:) = dcms.subj(s).dcmmat';
    end
    
    % Load all models into memory
    GCM = spm_dcm_load(P);    
    
elseif isfield(dcms,'gcmmat')
    GCM = load(dcms.gcmmat{1});
    GCM = GCM.GCM;
    ns = size(GCM,1);
    nm = size(GCM,2);
end

% Set timeseries or CSD estimation
for s = 1:ns
    for m = 1:nm
        if strcmpi(job.analysis,'CSD')
            GCM{s,m}.options.analysis = 'CSD';
        else
            if isfield(GCM{s,m},'options') && isfield(GCM{s,m},'analysis')
                GCM{s,m}.options = rmfield(GCM{s,m}.options,'analysis');
            end
        end
    end
end
    
% Estimate models if requested
switch job.est_type
    case EST_FULL_BMR
        GCM(:,1) = spm_dcm_fit(GCM(:,1));
        
        if nm > 1
            GCM = spm_dcm_bmr(GCM);
        end
    case EST_FULL_BMR_PEB
        GCM = spm_dcm_peb_fit(GCM);
    case EST_FULL
        GCM = spm_dcm_fit(GCM);    
    case EST_NONE
        % Do nothing
end

% Save
if isfield(job.output_type,'output_single')
    % Create single mat file
    dir  = job.output_type.output_single.dir{1};
    name = ['GCM_' job.output_type.output_single.name '.mat'];
    
    filename = fullfile(dir,name);
    save(filename,'GCM');
    
    out.gcmmat = {filename};
elseif (job.est_type ~= EST_NONE)
    % Update existing mat files
    for s = 1:ns
        for m = 1:nm
            DCM = GCM{s,m};
            save(P{s,m}, 'DCM');
        end
    end
    out.dcmmat = P;
end


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

%==========================================================================
function dep = vout_gcm_fmri(job)
%==========================================================================
if isfield(job.output_type,'output_single')
    dep(1)            = cfg_dep;
    dep(1).sname      = 'GCM mat File(s)';
    dep(1).src_output = substruct('.','gcmmat');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
else
    dep = [];
end