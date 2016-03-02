function estimate = spm_cfg_dcm_est
% SPM Configuration file for DCM estimation
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: spm_cfg_dcm_est.m 6735 2016-03-02 15:40:47Z peter $

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
subjects.name   = 'Per subject';
subjects.values = {subj};
subjects.help   = {'Create the subjects and select the models for each'};
subjects.num    = [1 Inf];

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
% output_single Output one .mat file for the group
%--------------------------------------------------------------------------
output_single         = cfg_branch;
output_single.tag     = 'single';
output_single.name    = 'Create group GCM_*.mat file';
output_single.val     = { dir name };
output_single.help    = {['Creates a single group-level DCM file ' ...
                          'containing a subjects x models cell array.']};
                      
% -------------------------------------------------------------------------
% output_overwrite_gcm Output a GCM file with existing name
%--------------------------------------------------------------------------
output_overwrite_gcm         = cfg_branch;
output_overwrite_gcm.tag     = 'overwrite_gcm';
output_overwrite_gcm.name    = 'Overwrite existing GCM_*.mat file';
output_overwrite_gcm.val     = {};
output_overwrite_gcm.help    = {['Overwrites existing group-level DCM file ' ...
                                 'with estimated models.']};

% -------------------------------------------------------------------------
% output_separate Output one .mat file per model
%--------------------------------------------------------------------------
output_separate         = cfg_branch;
output_separate.tag     = 'separate';
output_separate.name    = 'Overwrite existing individual DCM files';
output_separate.val     = {};
output_separate.help    = {'Updated existing individual DCM.mat files'};

% -------------------------------------------------------------------------
% output Choice of how many DCM.mat files to output
%--------------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output';
output.values  = { output_single output_overwrite_gcm output_separate };
output.val     = { output_single };
output.help    = {['Whether to create a single DCM file across all ' ...
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
dcms.val    = {gcmmat};
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
% analysis Analysis
% -------------------------------------------------------------------------
fmri_analysis         = cfg_menu;
fmri_analysis.tag     = 'analysis';
fmri_analysis.name    = 'Analysis';
fmri_analysis.labels  = {'default (time series)','cross-spectral densities'};
fmri_analysis.values  = {'time','csd'};
fmri_analysis.val     = {'time'};
fmri_analysis.help    = {['Whether to analyse in the time domain (for task-' ...
                     'based studies or stochastic DCM) or in the frequency '...
                     'domain (for resting state analysis with DCM for CSD']};
                 
fmri         = cfg_branch;
fmri.tag     = 'fmri';
fmri.name    = 'MRI specific options';
fmri.val     = {fmri_analysis};
                 
% -------------------------------------------------------------------------
% estimate Estimate
% -------------------------------------------------------------------------
estimate      = cfg_exbranch;
estimate.tag  = 'estimate';
estimate.name = 'DCM estimation';
estimate.val  = { dcms output est_type fmri };
estimate.help = {['Estimate the parameters and free energy (log model ' ...
                  'evidence) of first level DCMs for fMRI. Models ' ...
                  'are assembled into a Subjects x Models array and ' ...
                  'saved in group GCM_*.mat file']};
estimate.prog = @spm_run_dcm_est;
estimate.vout = @vout_dcm;

% -------------------------------------------------------------------------
% fmri Dynamic Causal Model for fMRI
% -------------------------------------------------------------------------
est         = cfg_choice; 
est.tag     = 'est';
est.name    = 'DCM estimation';
est.help    = {'Estimation of DCM models'};
est.values  = { estimate };

%==========================================================================
function out = spm_run_dcm_est(job)
%==========================================================================

dcms = job.dcms;

% Get selected estimation option
EST_FULL_BMR     = 1;
EST_FULL_BMR_PEB = 2;
EST_FULL         = 3;
EST_NONE         = 4;

est_type = job.est_type;

% Get selected input option
INPUT_DCM_BY_MODEL   = 1;
INPUT_DCM_BY_SUBJECT = 2;
INPUT_GCM            = 3;

if isfield(dcms,'model')
    input_type = INPUT_DCM_BY_MODEL;
elseif isfield(dcms,'subj')
    input_type = INPUT_DCM_BY_SUBJECT;
elseif isfield(dcms,'gcmmat')
    input_type = INPUT_GCM;
else
    error('Unknown input type');
end

% Get selected output option
OUTPUT_GCM_NEW       = 1;
OUTPUT_GCM_OVERWRITE = 2;
OUTPUT_DCM           = 3;

if isfield(job.output,'single')
    output_type = OUTPUT_GCM_NEW;
elseif isfield(job.output,'overwrite_gcm')
    output_type = OUTPUT_GCM_OVERWRITE;
elseif isfield(job.output,'separate')
    output_type = OUTPUT_DCM;
else
    error('Unknown output type');
end

% Validate options
if (input_type == INPUT_GCM && output_type == OUTPUT_DCM)
    error('If the input is a single GCM file, the output must be too');
end

if (output_type == OUTPUT_GCM_OVERWRITE && input_type ~= INPUT_GCM)
    error('To overwrite an existing GCM file, the input must be a GCM');
end

% Build subjects x models filename matrix
switch input_type
    case INPUT_DCM_BY_MODEL
        ns = length(dcms.model(1).dcmmat);
        nm = length(dcms.model);
        P  = cell(ns,nm);

        for m = 1:nm
            if length(dcms.model(m).dcmmat) ~= ns
                error(['Please ensure all models have the same number of ' ... 
                       'subjects']);
            end

            P(:,m) = dcms.model(m).dcmmat;
        end

        % Load all models into memory
        GCM = spm_dcm_load(P);    
    
    case INPUT_DCM_BY_SUBJECT
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
    
    case INPUT_GCM
        GCM = load(dcms.gcmmat{1});
        GCM = GCM.GCM;
        ns = size(GCM,1);
        nm = size(GCM,2);
end

% Set timeseries or CSD estimation (fMRI)
for s = 1:ns
    for m = 1:nm
        if strcmpi(job.fmri.analysis,'CSD')
            GCM{s,m}.options.analysis = 'CSD';
        else
            if isfield(GCM{s,m},'options') && isfield(GCM{s,m},'analysis')
                GCM{s,m}.options = rmfield(GCM{s,m}.options,'analysis');
            end
        end
    end
end

% Estimate models if requested
switch est_type
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
if output_type == OUTPUT_GCM_NEW
    % Create single mat file
    dir  = job.output.single.dir{1};
    name = ['GCM_' job.output.single.name '.mat'];
    
    filename = fullfile(dir,name);
    save(filename,'GCM');
    
    out.gcmmat = {filename};
elseif output_type == OUTPUT_GCM_OVERWRITE
    % Update existing gcm file
    filename = dcms.gcmmat{1};
    
    save(filename,'GCM');
    
    out.gcmmat = {filename};
else
    % Update existing mat files
    if (est_type ~= EST_NONE)
        for s = 1:ns
            for m = 1:nm
                DCM = GCM{s,m};
                F   = DCM.F;
                Ep  = DCM.Ep;
                Cp  = DCM.Cp;
                save(P{s,m}, 'DCM', 'F', 'Ep', 'Cp');
            end
        end
    end
    out.dcmmat = P;
end

%==========================================================================
function dep = vout_dcm(job)
%==========================================================================
if isfield(job.output,'single') || ...
        isfield(job.output,'overwrite_gcm')
    dep(1)            = cfg_dep;
    dep(1).sname      = 'GCM mat File(s)';
    dep(1).src_output = substruct('.','gcmmat');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
else
    dep(1)            = cfg_dep;
    dep(1).sname      = 'DCM mat File(s)';
    dep(1).src_output = substruct('.','dcmmat');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end