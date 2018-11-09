function fmri = spm_cfg_dcm_fmri
% SPM Configuration file for DCM for fMRI
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: spm_cfg_dcm_fmri.m 7479 2018-11-09 14:17:33Z peter $


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
% dcmmat_full Select template full model
% -------------------------------------------------------------------------
dcmmat_full      = dcmmat;
dcmmat_full.tag  = 'fulldcm';
dcmmat_full.name = 'Full DCM';
dcmmat_full.help = {['Select a DCM .mat file, specified using the GUI ' ...
    'for a single subject. This should be a ''full'' model with all ' ...
    'parameters of interest (e.g. connections) enabled.']};
dcmmat_full.num     = [1 1];

% -------------------------------------------------------------------------
% dcmmat_alt Select alternative template models
% -------------------------------------------------------------------------
dcmmat_alt = dcmmat;
dcmmat_alt.tag  = 'altdcm';
dcmmat_alt.name = 'Alternative DCMs';
dcmmat_alt.help = {['Select alternative DCMs, specified using the GUI ' ...
    'for a single subject. These should be reduced or nested versions of the ' ...
    'full template DCM, e.g. with some connections switched off (fixed ' ...
    'at their priors).']};
dcmmat_alt.num = [0 Inf];
dcmmat_alt.val = {''};
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
% region Select VOI_*.mat (for a single region across subjects)
% -------------------------------------------------------------------------
region         = cfg_files;
region.tag     = 'region';
region.name    = 'Region (VOI files)';
region.help    = {'Select all subjects'' VOI_*.mat for this region.'};
region.filter  = 'mat';
region.ufilter = '^VOI_.*\.mat$';
region.num     = [1 Inf];

% -------------------------------------------------------------------------
% regions Create set of regions
%--------------------------------------------------------------------------
vois       = cfg_repeat;
vois.tag    = 'vois';
vois.name   = 'Regions of interest';
vois.values = {region};
vois.help   = {'Select the regions of interest (VOIs) for all subjects'};
vois.num    = [1 Inf];

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
% spmmat Select SPM.mat for each subject
% -------------------------------------------------------------------------
spmmats      = spmmat;
spmmats.tag  = 'spmmats';
spmmats.name = 'SPM.mat files';
spmmats.help = {'Select all subjects'' SPM.mat files'};
spmmats.num  = [1 Inf];

% -------------------------------------------------------------------------
% session Session index
% -------------------------------------------------------------------------
session         = cfg_entry;
session.tag     = 'session';
session.name    = 'Which session';
session.help    = {'Enter the session (run) number.'};
session.strtype = 'n';
session.num     = [1 1];
session.val     = {1};

%--------------------------------------------------------------------------
% dir Directory
%--------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Output directory';
dir.help    = {'Select the directory where the output will be written.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% -------------------------------------------------------------------------
% name Model name
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output name';
name.help    = {['Specify a name for the group DCM file. The prefix ''GCM_'' ' ...
                'and suffix ''.mat'' are added automatically.']};
name.strtype = 's';
name.num     = [0 Inf];

% -------------------------------------------------------------------------
% val Val
% -------------------------------------------------------------------------
val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Condition';
val.help    = {['Values to include for this condition. Enter ''1'' ' ...
               'to include this condition (with no parameteric regressor). '...
               'Entering [1 0 1] would include this condition and '...
               'its second parametric regressor.']};
val.strtype = 'w';
val.num     = [1 Inf];

% -------------------------------------------------------------------------
% inp Inputs
% -------------------------------------------------------------------------
inp         = cfg_repeat;
inp.tag     = 'inputs';
inp.name    = 'Inputs';
inp.help    = {['Conditions to include and their parametric modulations (PMs). '...
               'Click ''New: Condition'' for each condition in your '...
               'SPM (i.e. SPM.U), up to the last condition you wish  '...
               'to include. For each Condition, enter:'] ...
               '1 (to include the condition in the DCM) or' ...
               '[1 1] (the condition and its 1st parameteric regressor) or' ...
               '[1 0 1] (the condition and its 2nd parameteric regressor)' ...
               'etc'};
inp.values  = { val };
inp.num     = [1 Inf];

% -------------------------------------------------------------------------
% subj Create single model
%--------------------------------------------------------------------------
model      = cfg_branch;
model.tag  = 'model';
model.name = 'Model';
model.val  = {dcmmat};
model.help = {'Corresponding model for each subject.'};

% -------------------------------------------------------------------------
% subjects Create set of models
%--------------------------------------------------------------------------
models        = cfg_repeat;
models.tag    = 'models';
models.name   = 'Per model';
models.values = {model};
models.help   = {'Select DCM.mat files per model.'};
models.num    = [1 Inf];

% -------------------------------------------------------------------------
% output Output branch (specify group batch)
% -------------------------------------------------------------------------
output         = cfg_branch;
output.tag     = 'output';
output.name    = 'Output';
output.val     = {dir name};
output.help    = {['Select where the group DCM (GCM) file containg all ' ...
                    'the DCMs'' filenames will be stored.']};

% -------------------------------------------------------------------------
% templates Templates branch (specify group batch)
% -------------------------------------------------------------------------
template       = cfg_branch;
template.tag   = 'template';
template.name  = 'Template DCMs';
template.val   = {dcmmat_full dcmmat_alt};
template.help  = {'Template DCMs to replicate over subjects.'};

% -------------------------------------------------------------------------
% data Design & data branch (specify group batch)
% -------------------------------------------------------------------------
data       = cfg_branch;
data.tag   = 'data';
data.name  = 'Design & data';
data.val   = {spmmats session vois};
data.help  = {'Experimental timing and timeseries for the DCMs.'};

% -------------------------------------------------------------------------
% group Specify group model space
% -------------------------------------------------------------------------
group      = cfg_exbranch;
group.tag  = 'group';
group.name = 'Specify group';
group.val  = {output template data};
group.help = {'Duplicates template DCM(s) for each subject.'...
    ['Before running this, create one or more example DCMs for one ' ...
    'subject using the Dynamic Causal Modelling button in the main ' ...
    'SPM window. Then use this batch to duplicate the DCM(s) for each ' ...
    'subject. The output is a file containing a cell array, with one ' ...
    'row per subject and one column per DCM (named GCM_<name>.mat).']};
group.prog = @spm_run_create_gcm;
group.vout = @vout_gcm_fmri;

% -------------------------------------------------------------------------
% regions Specify regions
% -------------------------------------------------------------------------
regions      = cfg_exbranch;
regions.tag  = 'regions';
regions.name = 'Region specification';
regions.val  = { dcmmat voimat };
regions.help = {'Insert new regions into a DCM model.'...
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
inputs.help = {'Insert new inputs into a DCM model.'...
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
fmri.values  = { group regions inputs };

%==========================================================================
function out = spm_run_create_gcm(job)
% Replicates template DCM(s) across subjects and creates a GCM file

% Unpack
out_dir  = job.output.dir;
out_name = job.output.name;
full_dcm = job.template.fulldcm;
alt_dcm  = job.template.altdcm;
spms     = job.data.spmmats;
sess     = job.data.session;
regions  = job.data.region;

nr = length(regions);     % number of regions
ns = length(spms);        % number of subjects
nm = length(alt_dcm) + 1; % number of models

% Validate length of VOI arrays
for r = 1:nr
    if length(regions{r}) ~= ns
        error('Pleasure ensure region %d has 1 VOI per subject',r);
    end
end

% Get template DCMs
templates    = cell(1,nm);
templates{1} = full_dcm{1};
for m = 1:length(alt_dcm)
    templates{m + 1} = alt_dcm{m};
end

% Get indices in the SPM of experimental inputs
DCM = load(templates{1});
DCM = DCM.DCM;
if isfield(DCM.U,'idx')
    idx = DCM.U.idx;
elseif isfield(DCM.U,'u') && size(DCM.U.u,2)==1 && ...
        strcmp(DCM.U.name{1},'null') ...
    % Resting state (backward compatibility)
    idx = 0;
else
    % GUI (backward compatibility)
    idx = get_conditions_from_ui(spms{1},sess); 
end

% Convert inputs to cell format for spm_dcm_U
inputs = {};
if ~isempty(idx) && (idx(1) ~= 0)
    max_cond = max(idx(:,1));
    inputs   = cell(1, max_cond);
    for i = 1:size(idx,1)        
        cond = idx(i,1);
        col  = idx(i,2);
        inputs{cond}(col) = 1;
    end   
end

GCM = cell(ns,nm);

for s = 1:ns
    % Get subject SPM dir
    subject_spm     = spms{s};      
    subject_spm_dir = fileparts(subject_spm);
    
    % Get subject VOIs
    subject_vois = cell(nr,1);
    for r = 1:nr
        subject_vois{r} = regions{r}{s};
    end
    
    for m = 1:nm
        % Copy template to subject
        new_dcm = fullfile(subject_spm_dir, ...
            sprintf('DCM_%s_m%04d.mat',out_name,m));
        copyfile(templates{m}, new_dcm);
        GCM{s,m} = new_dcm;
        
        % Update regions
        if ~isempty(inputs)
            spm_dcm_U(new_dcm,subject_spm,sess,inputs);
        end
        
        % Update timeseries
        spm_dcm_voi(new_dcm,subject_vois);
    end    
    
end

% Save GCM
out_file = fullfile(out_dir{1}, ['GCM_' out_name '.mat']);
save(out_file,'GCM');

out.gcmmat = {out_file};

%==========================================================================
function idx = get_conditions_from_ui(spmmat,sess)
% Prompts for condition names. Provides backward compatibility for DCMs
% specified prior to the introduction of DCM.U.idx 
% 
% idx - [n x 2] matrix for n conditions. The first column is the index in
%       SPM.U and the second is the regressor within the condition.

SPM = load(spmmat);
SPM = SPM.SPM;

Sess = SPM.Sess(sess);
u    = length(Sess.U);

if isempty(Sess.U)
    idx = 0;
    return;
end

idx  = [];
for i = 1:u
    for j = 1:length(Sess.U(i).name)
        str = ['include ' Sess.U(i).name{j} '?'];
        if spm_input(str,'+1','y/n',[1 0],1)
            idx = [idx; i j];
        end
    end
end

if isempty(idx), idx = 0; end

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
function dep = vout_gcm_fmri(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'GCM mat File(s)';
dep(1).src_output = substruct('.','gcmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
