function second_level = spm_cfg_dcm_peb
% SPM Configuration file for second-level DCM (PEB)
%__________________________________________________________________________
% Copyright (C) 2016-2017 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_cfg_dcm_peb.m 7007 2017-02-07 10:15:24Z guillaume $


%==========================================================================
% Directory / filename selection
%==========================================================================

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

%--------------------------------------------------------------------------
% name Model name
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Specify a name for the output.'};
name.strtype = 's';
name.num     = [0 Inf];

%==========================================================================
% DCM / PEB file selection
%==========================================================================

%--------------------------------------------------------------------------
% dcmmat Select DCM_*.mat
%--------------------------------------------------------------------------
dcmmat         = cfg_files;
dcmmat.tag     = 'dcmmat';
dcmmat.name    = 'Select DCM files';
dcmmat.help    = {'Select DCM_*.mat files.'};
dcmmat.filter  = 'mat';
dcmmat.ufilter = '^DCM_.*\.mat$';
dcmmat.num     = [1 Inf];

%--------------------------------------------------------------------------
% subj Create single subject
%--------------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {dcmmat};
subj.help = {'Subject with one or more models.'};
    
%--------------------------------------------------------------------------
% peb_mat Select PEB_*.mat
%--------------------------------------------------------------------------
peb_mat         = cfg_files;
peb_mat.tag     = 'peb_mat';
peb_mat.name    = 'Select PEB file';
peb_mat.help    = {'Select PEB_*.mat file.'};
peb_mat.filter  = 'mat';
peb_mat.ufilter = '^PEB_.*\.mat$';
peb_mat.num     = [1 1];

%--------------------------------------------------------------------------
% peb_mat Select model_space_*.mat
%--------------------------------------------------------------------------
model_space_mat         = cfg_files;
model_space_mat.tag     = 'model_space_mat';
model_space_mat.name    = 'DCMs';
model_space_mat.help    = {['Select group DCM file (GCM_*.mat). This is a ' ...
                           'cell array with one row per subject and one ' ...
                           'column per DCM.']};
model_space_mat.filter  = 'mat';
model_space_mat.ufilter = '^GCM.*\.mat$';
model_space_mat.num     = [1 Inf];

%==========================================================================
% Covariates entry
%==========================================================================

%--------------------------------------------------------------------------
% design Design matrix
%--------------------------------------------------------------------------
cov_design         = cfg_entry;
cov_design.tag     = 'cov_design';
cov_design.name    = 'Design matrix';
cov_design.help    = {['Enter or paste the N x C design matrix for N ' ...
                      'subjects and C covariates. Note that a column of '...
                      'ones will automatically be added to the start of the '...
                      'matrix, to model the group mean.']};
cov_design.strtype = 'r';
cov_design.num     = [Inf Inf];

%--------------------------------------------------------------------------
% cov_name Name for a column in the design matrix
%--------------------------------------------------------------------------
cov_name         = cfg_entry;
cov_name.tag     = 'name';
cov_name.name    = 'Name';
cov_name.help    = {'Enter a name for a covariate.'};
cov_name.strtype = 's';
cov_name.num     = [0 Inf];

%--------------------------------------------------------------------------
% cov_names Contains the names for the covariates
%--------------------------------------------------------------------------
cov_names         = cfg_repeat;
cov_names.tag     = 'names';
cov_names.name    = 'Covariate names';
cov_names.values  = { cov_name };
cov_names.num     = [0 Inf];
cov_names.help    = {['Enter names for each covariate (excluding the mean ' ...
                     'regressor which is added automatically).']};

%--------------------------------------------------------------------------
% design_mtx Specify whole design matrix 
%--------------------------------------------------------------------------
design_mtx         = cfg_branch;
design_mtx.tag     = 'design_mtx';
design_mtx.name    = 'Specify design matrix';
design_mtx.val     = { cov_design cov_names };
design_mtx.help    = {'Specify the second-level design matrix.'};

%--------------------------------------------------------------------------
% cov_val Value
%--------------------------------------------------------------------------
cov_val         = cfg_entry;
cov_val.tag     = 'value';
cov_val.name    = 'Value';
cov_val.help    = {'Enter the vector of regressor values, one element ' ...
                   'per subject.'};
cov_val.strtype = 'r';
cov_val.num     = [Inf 1];

%--------------------------------------------------------------------------
% covariate A single covariate
%--------------------------------------------------------------------------
regressor         = cfg_branch;
regressor.tag     = 'regressor';
regressor.name    = 'Covariate';
regressor.val     = {cov_name cov_val };
regressor.help    = {'Regressor.'};

%--------------------------------------------------------------------------
% covariate Specify design matrix per covariate
%--------------------------------------------------------------------------
regressors         = cfg_repeat;
regressors.tag     = 'regressors';
regressors.name    = 'Specify covariates individually';
regressors.values  = { regressor };
regressors.help    = {'Specify the second-level design matrix one '...
                     'covariate (regressor) at a time. Note that a ' ...
                     'column of ones to model the mean across subjects ' ...
                     'is added automatically.'};
regressors.num     = [1 Inf];

%--------------------------------------------------------------------------
% none No covariates
%--------------------------------------------------------------------------
cov_none         = cfg_branch;
cov_none.tag     = 'none';
cov_none.name    = 'None';
cov_none.val     = {};
cov_none.help    = {'Include no covariates (only the group mean for each '...
                'connection).'};

%--------------------------------------------------------------------------
% covariates Covariates branch
%--------------------------------------------------------------------------
covariates         = cfg_choice;
covariates.tag     = 'cov';
covariates.name    = 'Covariates';
covariates.values  = { cov_none design_mtx regressors };
covariates.help    = {['Specify between-subjects effects (covariates). The ' ...
                      'covariates may be entered all at once as a design ' ...
                      'matrix or individually. Note that if performing ' ...
                      'bayesian model comparison, only the first covariate ' ...
                      'will be treated as being of experimental interest.'] '' ...                      
                      ['Each parameter in the estimated PEB model will ' ...
                      'represent the influence of a covariate on a ' ...
                      'connection. If none is set, only the group mean will '...
                      'be estimated.']};
covariates.val    = {cov_none};
    
%==========================================================================
% PEB fields selection
%==========================================================================

%--------------------------------------------------------------------------
% field_default Select fields A,B
%--------------------------------------------------------------------------
field_default      = cfg_const;
field_default.tag  = 'default';
field_default.name = 'A- and B-matrix';
field_default.val  = {{'A','B'}};
field_default.help = {'A- and B-matrix.'};

%--------------------------------------------------------------------------
% field_all Select all fields
%--------------------------------------------------------------------------
field_all      = cfg_const;
field_all.tag  = 'all';
field_all.name = 'All';
field_all.val  = {'All fields'};
field_all.help = {'All fields.'};

%--------------------------------------------------------------------------
% field_entry Custom field entry
%--------------------------------------------------------------------------
field_entry  = cfg_entry;
field_entry.name = 'Enter manually';
field_entry.tag  = 'custom';
field_entry.help = {['Enter the fields as a cell array e.g. {''A''} ' ...
                     'or {''A'', ''C''} or {''B(:,:,1)'', ''B(:,:,3)''}']};
field_entry.strtype = 'e';
field_entry.num     = [0 Inf];

%--------------------------------------------------------------------------
% fields DCM fields to include
%--------------------------------------------------------------------------
fields        = cfg_choice;
fields.tag    = 'fields';
fields.name   = 'Fields';
fields.values = {field_default field_all field_entry};
fields.help   = {'Select the fields of the DCM to include in the model.' '' ...
                  'A- and B-matrix: Includes all A- and B- connections' ...
                  'All: Includes all fields' ...
                  'Enter manually: Enter a cell array e.g. {''A'',''C''}'};
fields.val    = {field_default};

%==========================================================================
% DCM model index selection
%==========================================================================

%--------------------------------------------------------------------------
% dcm_all Select all DCMs
%--------------------------------------------------------------------------
dcm_all      = cfg_const;
dcm_all.tag  = 'all';
dcm_all.name = 'All';
dcm_all.val  = {'All DCMs'};
dcm_all.help = {'All DCMs.'};

%--------------------------------------------------------------------------
% dcm_idx Single DCM index selection
%--------------------------------------------------------------------------
dcm_sel_idx  = cfg_entry;
dcm_sel_idx.name = 'Selected DCM index';
dcm_sel_idx.tag  = 'index';
dcm_sel_idx.help = {['Select index of the DCM (within subject) on which to ' ...
                 'to build the PEB - e.g. 1 to use each subject''s ' ...
                 'first DCM.']};
dcm_sel_idx.strtype = 'r';
dcm_sel_idx.num     = [1 Inf];

%--------------------------------------------------------------------------
% dcms Which DCMs to include
%--------------------------------------------------------------------------
dcm_idx_1     = dcm_sel_idx;
dcm_idx_1.val = {1};

dcm_idx        = cfg_choice;
dcm_idx.tag    = 'dcm';
dcm_idx.name   = 'DCM index';
dcm_idx.values = {dcm_idx_1 dcm_all};
dcm_idx.help   = {['If each subject has multiple DCMs, select which DCM to use ' ...
                  'for this analysis (or select all). For a standard PEB ' ... 
                  'analysis, it is recommended to specify the PEB over the ' ...
                  'first DCM only, then compare models at the second level.']};
dcm_idx.val    = {dcm_idx_1};
           
%==========================================================================
% Priors on log precision (between-subjects variability) entry
%==========================================================================

%--------------------------------------------------------------------------
% priors_log_precision_mu Priors on log precision expectation
%--------------------------------------------------------------------------
priors_log_precision_mu      = cfg_entry;
priors_log_precision_mu.name = 'Expectation';
priors_log_precision_mu.tag  = 'expectation';
priors_log_precision_mu.help = {['Prior expectation of the log precision ' ...
      'parameters (M.hE), which scale each precision component. The default ' ...
      'is 0, which is translated internally to exp(0) = 1.']};
priors_log_precision_mu.strtype = 'r';
priors_log_precision_mu.num     = [1 1];
priors_log_precision_mu.val     = {0};

%--------------------------------------------------------------------------
% priors_log_precision_var Priors on log precision variance
%--------------------------------------------------------------------------
priors_log_precision_var      = cfg_entry;
priors_log_precision_var.name = 'Uncertainty';
priors_log_precision_var.tag  = 'var';
priors_log_precision_var.help = {['Uncertainty over the prior expectation ' ...
    'of the log precision parameters (M.hC), which scale each precision ' ...
    'component. The default is 1/16.']};
priors_log_precision_var.strtype = 'r';
priors_log_precision_var.num     = [1 1];
priors_log_precision_var.val     = {1/16};

%--------------------------------------------------------------------------
% group priors_parameters_ratio Priors on log precision variance
%--------------------------------------------------------------------------
group_priors_parameters_ratio      = cfg_entry;
group_priors_parameters_ratio.name = 'Group ratio';
group_priors_parameters_ratio.tag  = 'group_ratio';
group_priors_parameters_ratio.help = {['The group ratio (M.alpha) expresses ' ...
     'our prior uncertainty about second-level (group ' ...
     'level GLM) parameters. The default is 1, meaning our uncertainty about ' ...
     'the connection strengths at the second level is the same as our ' ...
     'uncertainty about the connection strengths at the first level. '] '' ...
     ['Internally, the prior covariance of second level parameters is set ' ...
     'to DCM{1}.M.pC / alpha, where DCM{1} is the first DCM provided.']};
group_priors_parameters_ratio.strtype = 'r';
group_priors_parameters_ratio.num     = [1 1];
group_priors_parameters_ratio.val     = {1};

%--------------------------------------------------------------------------
% priors_parameters_ratio Priors on log precision variance
%--------------------------------------------------------------------------
priors_parameters_ratio  = cfg_entry;
priors_parameters_ratio.name = 'Within:between ratio';
priors_parameters_ratio.tag  = 'ratio';
priors_parameters_ratio.help = {'Within:between variance ratio (M.beta). ' ...
    ['This ratio controls the expected between-subjects ' ...
    'variability for each connection (DCM parameter). The default is 16, ' ...
    'meaning we expect the variability in connection strengths across ' ...
    'subjects to be 1/16 of our uncertainty about connection strengths at ' ...
    'the first level.'] '' ...
    ['Internally, the prior second level covariance of DCM parameters (M.pC) ' ...
    'is set to DCM{1}.M.pC / M.beta, where DCM{1} is the first DCM provided.']};
priors_parameters_ratio.strtype = 'r';
priors_parameters_ratio.num     = [1 1];
priors_parameters_ratio.val     = {16};


%--------------------------------------------------------------------------
% priors_between Priors on log precision branch
%--------------------------------------------------------------------------
priors_between      = cfg_branch;
priors_between.tag  = 'priors_between';
priors_between.name = 'Between-subjects variability';
priors_between.val  = { priors_parameters_ratio ...
                        priors_log_precision_mu ...
                        priors_log_precision_var};
priors_between.help = {['Between-subjects variability over second-' ...
     'level parameters.'], '' ...
     ['A multi-component model is used. Each component is a ' ...
      '[p x p] precision matrix given p DCM parameters, where elements on ' ...
      'the diagonal represent the precision (inverse variance) across ' ...
      'subjects of each DCM connection. These precisions are set via the ' ...
      'Within:between ratio, below. Each precision component is scaled by ' ...
      'a hyper-parameter, which is estimated from the data. The prior ' ...
      'expectation and uncertainty of these hyper-parameters are also '...
      'set below.']};
  

%--------------------------------------------------------------------------
% priors_third Priors on log precision branch
%--------------------------------------------------------------------------
priors_glm         = cfg_branch;
priors_glm.tag     = 'priors_glm';
priors_glm.name    = 'Second level (GLM) priors';
priors_glm.val     = { group_priors_parameters_ratio};
priors_glm.help    = {['Priors on the expected values of the second ' ...
                       'level General Linear Model (GLM) parameters.']};

%--------------------------------------------------------------------------
% show_review Select whether to review results
%--------------------------------------------------------------------------
show_review        = cfg_menu;
show_review.tag    = 'show_review';
show_review.name   = 'Review PEB parameters';
show_review.labels = {'Yes','No'};
show_review.values = {1,0};
show_review.val    = {1};
show_review.help   = {'Review PEB parameters'};

%==========================================================================
% PEB specification batch
%==========================================================================
% Set show review default to off
sr      = show_review;
sr.val  = {0};

specify      = cfg_exbranch;
specify.tag  = 'specify';
specify.name = 'Specify / Estimate PEB';
specify.val  = { name model_space_mat dcm_idx covariates fields ...
                 priors_between priors_glm sr };
specify.help = {['Specifies and estimates a second-level DCM (PEB) model. ' ...
                 'A PEB model will be created for each first level DCM.' ]};
            
specify.prog = @spm_run_create_peb;
specify.vout = @vout_peb;

%==========================================================================
% PEB reduce / average / compare batch
%==========================================================================
model_space_mat_op = model_space_mat;
model_space_mat_op.num = [0 Inf];
model_space_mat_op.val = {''};

null_prior_covariance  = cfg_entry;
null_prior_covariance.name = 'Null prior covariance';
null_prior_covariance.tag  = 'nullpcov';
null_prior_covariance.help = {'gamma'};
null_prior_covariance.strtype = 'r';
null_prior_covariance.num     = [1 1];
null_prior_covariance.val     = {1/16};

peb_compare      = cfg_exbranch;
peb_compare.tag  = 'compare';
peb_compare.name = 'Compare / Average PEB models';
peb_compare.val  = { peb_mat model_space_mat_op null_prior_covariance show_review};
peb_compare.help = {['Addresses the question: which combination of ' ...
    'connections best explains the commonalities across subjects and ' ...
    'the group differences between subjects?'] '' ...
    ['If a group difference is to be investigated, this should be the ' ...
     'first manually entered between-subjects covariate in the PEB.'] ''  ...
    ['This analysis is performed by comparing a PEB model to nested ' ...
    'sub-models where ' ...
    'certain parameters have been disabled (fixed at their prior mean of ' ...
    'zero). Parameters are then averaged over reduced models to give an '...
    'averaged PEB (referred to as a Bayesian Model Average, BMA). Each ' ...
    'parameter in the PEB or BMA represents the effect of one between-' ...
    'subjects covariate on one connection.'] '' ...
    ['If only one first-level DCM is provided per subject, a search is made ' ...
    'over reduced PEB models to prune away any parameters not contributing ' ... 
    'to the model evidence. If multiple DCMs are provided per subject, ' ...
    'these are used to define the combinations of second level parameters ' ...
    'to be compared.']};
peb_compare.prog = @spm_run_bmc;
peb_compare.vout = @vout_bma;

reduce_all       = cfg_exbranch;
reduce_all.tag   = 'reduce_all';
reduce_all.name  = 'Search nested PEB models';
reduce_all.val   = { peb_mat model_space_mat_op null_prior_covariance show_review};
reduce_all.help  = {['Optimises a PEB model by trying different ' ...
                         'combinations of switching off parameters (fixing ' ...
                         'them at their prior value), where doing so does ' ...
                         'not reduce the model evidence. Any parameters ' ...
                         'not contributing will be set to zero.'] '' ...
                         ['This is equivilant to the function of '... 
                         'spm_dcm_post_hoc but on the second level (PEB) ' ...
                         'parameters.']};
reduce_all.prog = @spm_run_bmr_all;
reduce_all.vout = @vout_bma;

%==========================================================================
% PEB review batch
%==========================================================================
review      = cfg_exbranch;
review.tag  = 'peb_review';
review.name = 'Review PEB';
review.val  = { peb_mat model_space_mat_op };
review.help = {'Reviews PEB results'};
review.prog = @spm_run_dcm_peb_review;

%==========================================================================
% PREDICT leave-one-out cross validation
%==========================================================================

covariates_min1 = covariates;
covariates_min1.val{1} = covariates_min1.values{2};
covariates_min1.values(1) = [];

predict      = cfg_exbranch;
predict.tag  = 'predict';
predict.name = 'Predict (cross-validation)';
predict.val  = { name model_space_mat dcm_idx covariates_min1 fields ...
                 priors_between priors_glm };
predict.help = {['Builds a PEB model on all but one subjects, and uses ' ...
                 'it to predict a between-subjects effect, such as ' ...
                 'group membership, in the remaining subject. This process ' ...
                 'is repeated, leaving out each subject in turn ' ...
                 '(leave-one-out cross-validation) to establish the cross- ' ...
                 'validation accuracy of the model.'] ...
                 ['The first covariate (after the automatically inserted ' ...
                  'mean regressor) is used as the predictor variable.']};            
predict.prog = @spm_run_dcm_loo;
predict.vout = @vout_loo;

%==========================================================================
% second_level Second level DCM batch
%==========================================================================
second_level         = cfg_choice; 
second_level.tag     = 'peb';
second_level.name    = 'Second level';
second_level.help    = {'Parametric Empirical Bayes for DCM.'};
second_level.values  = { specify peb_compare reduce_all review predict };


%==========================================================================
function dep = vout_bma(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'BMA mat filename';
dep(1).src_output = substruct('.','bmamat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

%==========================================================================
function dep = vout_loo(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'LOO mat File(s)';
dep(1).src_output = substruct('.','loo_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

%==========================================================================
function dep = vout_peb(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'PEB mat File(s)';
dep(1).src_output = substruct('.','peb_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});


%==========================================================================
function out = spm_run_dcm_peb_review(job)
%==========================================================================
% Run the PEB review batch
P   = job.peb_mat;
DCM = job.model_space_mat;
spm_dcm_peb_review(P{1},DCM);
out = job.peb_mat;

%==========================================================================
function out = spm_run_create_peb(job)
%==========================================================================
% Run the PEB specification / estimation batch

[GCM,M,field,gcm_file] = prepare_peb_inputs(job);
    
% Specify / estimate PEB on full model only
dir_out = fileparts(gcm_file);
name    = job.name;    
PEB     = spm_dcm_peb(GCM,M,field);

% Write PEB
peb_filename = fullfile(dir_out,['PEB_' name '.mat']);
save(peb_filename,'PEB', spm_get_defaults('mat.format'));

% Review PEB
if job.show_review == 1
    spm_dcm_peb_review(peb_filename,GCM);
end

out.peb_mat = {peb_filename};

%==========================================================================
function out = spm_run_dcm_loo(job)
%==========================================================================
% Run leave-one-out cross validation

[GCM,M,field,gcm_file] = prepare_peb_inputs(job);

[qE,qC,Q] = spm_dcm_loo(GCM,M,field);

% Write output
dir_out      = fileparts(gcm_file);
name         = job.name;    
loo_filename = fullfile(dir_out,['LOO_' name '.mat']);

save(loo_filename,'qE','qC','Q', spm_get_defaults('mat.format'));

out.loo_mat = {loo_filename};

%==========================================================================
function [GCM,M,field,gcm_file] = prepare_peb_inputs(job)
%==========================================================================
% Prepare the inputs needed to specify a PEB or run LOO

[GCM,gcm_file] = load_dcm(job);

ns = size(GCM,1);

if ~isfield(GCM{1},'Ep')
    error('Please estimate DCMs before second-level analysis.');
end

% DCM field(s)
if isfield(job.fields,'default')
    field = {'A','B'};
elseif isfield(job.fields,'all')
    field = 'all';
else
    field = job.fields.custom;
end

Xnames = {'Group mean'};

X = ones(ns,1);

% Covariates
if isfield(job.cov, 'none')
    % Do nothing
        
elseif isfield(job.cov, 'design_mtx')
    % Whole design matrix entered
    x = job.cov.design_mtx.cov_design;
    
    if size(x,1) ~= ns
        error('Please ensure design matrix has one row per subject.');
    end
    
    X = [X x];
    
    Xnames = [Xnames job.cov.design_mtx.name];
elseif isfield(job.cov, 'regressor')
    % Design matrix entered per-regressor
    
    regressors = job.cov.regressor;
       
    for r = 1:length(regressors)
        regressor = regressors(r).value;
        name      = regressors(r).name;
        
        if size(regressor,1) ~= ns
            error('Please ensure regressor %d has one row per subject.',r);
        end
        
        X = [X regressor];
        Xnames = [Xnames name];
    end
end

% Ensure a mean column wasn't entered accidently
if size(X,2) > 1
    bad = find(~any(diff(X(:,2:end))));
    if ~isempty(bad)
        error('Please check regressor %d.', bad);
    end
end

if size(X,2) ~= length(Xnames)
    error('Please ensure there is one covariate name per covariate.');
end

% Priors / covariance components
M = struct();
M.alpha  = job.priors_glm.group_ratio;
M.beta   = job.priors_between.ratio;
M.hE     = job.priors_between.expectation;
M.hC     = job.priors_between.var;
M.Q      = 'single';
M.X      = X;
M.Xnames = Xnames;

%==========================================================================
function out = spm_run_bmr_all(job)
%==========================================================================
% Run a search over all reduced PEB models and average

GCM = load_dcm(job);
nm  = size(GCM,2);

if nm ~= 1
    disp('Running search on the full DCM only.');        
end

GCM = GCM(:,1);

out = run_peb_bmc_internal(job,GCM);

%==========================================================================
function out = spm_run_bmc(job)
%==========================================================================
% Run a model comparison between specified reduced PEB models

GCM = load_dcm(job);
nm  = size(GCM,2);

if nm < 2
    error('More than one DCM per subject is required for BMC.');
end

out = run_peb_bmc_internal(job,GCM);

%==========================================================================
function out = run_peb_bmc_internal(job, GCM)
%==========================================================================
% Search / model comparison both use spm_dcm_peb_bmc, which is called here

PEB = load(job.peb_mat{1});
PEB = PEB.PEB;

PEB.gamma = job.nullpcov;

nm  = size(GCM,2);

% Run BMA on defined reduced models or all submodels
if nm > 1
    BMA = spm_dcm_peb_bmc(PEB,GCM(1,:));
else
    BMA = spm_dcm_peb_bmc(PEB);
end

% Write BMA
[dir_out, name] = fileparts(job.peb_mat{1});
filename = fullfile(dir_out, ['BMA_' name '.mat']);
save(filename,'BMA', spm_get_defaults('mat.format'));

out.bmamat = filename;

% Review BMA
if job.show_review == 1
    
    DCM = job.model_space_mat;
    if ~isempty(DCM)
        DCM = load(DCM{1});
        if isfield(DCM,'GCM')
            DCM = DCM.GCM;
        else
            DCM = DCM.DCM;
        end
    end    
    
    spm_dcm_peb_review(BMA,DCM);
end

%==========================================================================
function [GCM,gcm_file] = load_dcm(job)
%==========================================================================
% Load and validate selected model space

gcm_file = char(job.model_space_mat);
GCM      = load(gcm_file);
if ~isfield(GCM,'GCM')
    error('Provided file is not a valid model space.');
end
GCM = GCM.GCM;

% Limit to specific model(s) if requested
if isfield(job,'dcm') && isfield(job.dcm,'index')
    index = job.dcm.index;
    GCM = GCM(:,index);
end
