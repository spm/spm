function second_level = spm_cfg_dcm_peb
% SPM Configuration file for second-level DCM (PEB)
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_cfg_dcm_fmri.m 6565 2015-09-30 10:42:14Z peter $


% =========================================================================
% Directory / filename selection
% =========================================================================

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
name.help    = {'Specify a name for the output'};
name.strtype = 's';
name.num     = [0 Inf];

% =========================================================================
% DCM / PEB file selection
% =========================================================================

% -------------------------------------------------------------------------
% dcmmat Select DCM_*.mat
% -------------------------------------------------------------------------
dcmmat         = cfg_files;
dcmmat.tag     = 'dcmmat';
dcmmat.name    = 'Select DCM files';
dcmmat.help    = {'Select DCM_*.mat files.'};
dcmmat.filter  = 'mat';
dcmmat.ufilter = '^DCM_.*\.mat$';
dcmmat.num     = [1 Inf];

% -------------------------------------------------------------------------
% subj Create single subject
%--------------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {dcmmat};
subj.help = {'Subject with one or more models.'};
    
% -------------------------------------------------------------------------
% pebmat Select PEB_*.mat
% -------------------------------------------------------------------------
pebmat         = cfg_files;
pebmat.tag     = 'pebmat';
pebmat.name    = 'Select PEB file';
pebmat.help    = {'Select PEB_*.mat file to review.'};
pebmat.filter  = 'mat';
pebmat.ufilter = '^PEB_.*\.mat$';
pebmat.num     = [1 1];

% -------------------------------------------------------------------------
% pebmat Select model_space_*.mat
% -------------------------------------------------------------------------
model_space_mat         = cfg_files;
model_space_mat.tag     = 'model_space_mat';
model_space_mat.name    = 'First-level DCMs';
model_space_mat.help    = {'Select model_space_*.mat file or one DCM ' ...
                           'per subject'};
model_space_mat.filter  = 'mat';
model_space_mat.ufilter = '^GCM.*\.mat$';
model_space_mat.num     = [1 Inf];

% =========================================================================
% Covariates entry
% =========================================================================

% ---------------------------------------------------------------------
% design Design matrix
% ---------------------------------------------------------------------
cov_design         = cfg_entry;
cov_design.tag     = 'cov_design';
cov_design.name    = 'Design matrix';
cov_design.help    = {['Enter or paste the N x C design matrix for N ' ...
                      'subjects and C covariates. Note that a column of '...
                      'ones will automatically be added to the start of the '...
                      'matrix, to model the group mean.']};
cov_design.strtype = 'r';
cov_design.num     = [Inf Inf];

% ---------------------------------------------------------------------
% cov_name Name for a column in the design matrix
% ---------------------------------------------------------------------
cov_name         = cfg_entry;
cov_name.tag     = 'cov_name';
cov_name.name    = 'Name';
cov_name.help    = {'Enter a name for a covariate'};
cov_name.strtype = 's';
cov_name.num     = [0 Inf];

% ---------------------------------------------------------------------
% cov_names Contains the names for the covariates
% ---------------------------------------------------------------------
cov_names         = cfg_repeat;
cov_names.tag     = 'cov_names';
cov_names.name    = 'Covariate names';
cov_names.values  = { cov_name };
cov_names.num     = [0 Inf];
cov_names.help   = {['Enter names for each covariate (excluding the mean ' ...
                     'regressor which is added automatically).']};

% ---------------------------------------------------------------------
% design_mtx Specify whole design matrix 
% ---------------------------------------------------------------------
design_mtx         = cfg_branch;
design_mtx.tag     = 'design_mtx';
design_mtx.name    = 'Specify design matrix';
design_mtx.val     = { cov_design cov_names };
design_mtx.help    = {'Specify the second-level design matrix.'};

% ---------------------------------------------------------------------
% cov_val Value
% ---------------------------------------------------------------------
cov_val         = cfg_entry;
cov_val.tag     = 'cov_val';
cov_val.name    = 'Value';
cov_val.help    = {'Enter the vector of regressor values, one element ' ...
                   'per subject.'};
cov_val.strtype = 'r';
cov_val.num     = [Inf 1];

% ---------------------------------------------------------------------
% covariate A single covariate
% ---------------------------------------------------------------------
regressor         = cfg_branch;
regressor.tag     = 'regressor';
regressor.name    = 'Covariate';
regressor.val     = {cov_name cov_val };
regressor.help    = {'regressor'};

% ---------------------------------------------------------------------
% covariate Specify design matrix per covariate
% ---------------------------------------------------------------------
regressors         = cfg_repeat;
regressors.tag     = 'regressors';
regressors.name    = 'Specify covariates individually';
regressors.values  = { regressor };
regressors.help    = {'Specify the second-level design matrix one '...
                     'covariate (regressor) at a time. Note that a ' ...
                     'column of ones to model the mean across subjects ' ...
                     'is added automatically.'};
regressors.num     = [1 Inf];

% ---------------------------------------------------------------------
% none No covariates
% ---------------------------------------------------------------------
cov_none         = cfg_branch;
cov_none.tag     = 'none';
cov_none.name    = 'None';
cov_none.val     = {};
cov_none.help    = {'Include no covariates (only the group mean for each '...
                'connection)'};

% ---------------------------------------------------------------------
% covariates Covariates branch
% ---------------------------------------------------------------------
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
    
% =========================================================================
% PEB fields selection
% =========================================================================

% ---------------------------------------------------------------------
% field_default Select fields A,B
% ---------------------------------------------------------------------
field_default  = cfg_const;
field_default.tag  = 'field_default';
field_default.name = 'A- and B-matrix';
field_default.val = {{'A','B'}};

% ---------------------------------------------------------------------
% field_all Select all fields
% ---------------------------------------------------------------------
field_all  = cfg_const;
field_all.tag  = 'field_all';
field_all.name = 'All';
field_all.val = {'All fields'};

% ---------------------------------------------------------------------
% field_entry Custom field entry
% ---------------------------------------------------------------------
field_entry  = cfg_entry;
field_entry.name = 'Enter manually';
field_entry.tag  = 'field_entry';
field_entry.help = {'Enter the fields e.g. A or A,C'};
field_entry.strtype = 'e';
field_entry.num     = [0 Inf];

% ---------------------------------------------------------------------
% fields DCM fields to include
% ---------------------------------------------------------------------
fields         = cfg_choice;
fields.tag    = 'fields';
fields.name   = 'Fields';
fields.values = {field_default field_all field_entry};
fields.help   = {'Select the fields of the DCM to include in the model.' '' ...
                  'A- and B-matrix: Includes all A- and B- connections' ...
                  'All: Includes all fields' ...
                  'Enter manually: Enter a cell array e.g. {''A'',''C''}'};
fields.val    = {field_default};
                       
% =========================================================================
% Bayesian Model Comparison selection
% =========================================================================

bmc         = cfg_menu;
bmc.tag     = 'bmc';
bmc.name    = 'Bayesian Model Comparison';
bmc.labels  = {'Yes','No'};
bmc.values  = {1,0};
bmc.val     = {1};
bmc.help    = {['If set to Yes, a comparison is performed to identify which ' ...
                'first level DCM best explains between-subjects effects. '...
                'If set to No, a PEB model for each DCM is simply '...
                'returned'] '' ...
               ['The comparison works by estimating a PEB model for each ' ...
                'DCM and each combination of second level covariates. '...
                'The PEB with the greatest evidence over the joint ' ...
                'space of first and second-level models is returned.']};

            
% =========================================================================
% Priors on log precision (between-subjects variability) entry
% =========================================================================

% ---------------------------------------------------------------------
% priors_log_precision_mu Priors on log precision expectation
% ---------------------------------------------------------------------
priors_log_precision_mu  = cfg_entry;
priors_log_precision_mu.name = 'Expectation';
priors_log_precision_mu.tag  = 'expectation';
priors_log_precision_mu.help = {['Prior expectation of the log precision ' ...
      'parameters (M.hE), which scale each precision component. The default ' ...
      'is log(0) = 1.']};
priors_log_precision_mu.strtype = 'r';
priors_log_precision_mu.num     = [1 1];
priors_log_precision_mu.val     = {0};

% ---------------------------------------------------------------------
% priors_log_precision_var Priors on log precision variance
% ---------------------------------------------------------------------
priors_log_precision_var  = cfg_entry;
priors_log_precision_var.name = 'Uncertainty';
priors_log_precision_var.tag  = 'var';
priors_log_precision_var.help = {['Uncertainty over the prior expectation ' ...
    'of the log precision parameters (M.hC), which scale each precision ' ...
    'component. The default is 1/16.']};
priors_log_precision_var.strtype = 'r';
priors_log_precision_var.num     = [1 1];
priors_log_precision_var.val     = {1/16};

% ---------------------------------------------------------------------
% components Covariance component selection
% ---------------------------------------------------------------------
components = cfg_menu;
components.tag    = 'components';
components.name   = 'Precision components';
components.labels = {'Single','Fields','All','None'};
components.values = {'single','fields','all','none'};
components.val    = {'single'};
components.help   = {'Select the precision components to create:' ...
    ['Single - One precision component is specified, ' ...
    'assuming all parameters vary equally over subjects.'] ...
    ['Fields - One precision component is specified for each field ' ...
     '(A,B etc).'] ...
    'All - One precision component is specified for each DCM connection.' ...
    ['None - No precision components are specified, thus all parameters ' ...
     'are fixed effects.']};

% ---------------------------------------------------------------------
% priors_parameters_ratio Priors on log precision variance
% ---------------------------------------------------------------------
priors_parameters_ratio  = cfg_entry;
priors_parameters_ratio.name = 'Within:between ratio';
priors_parameters_ratio.tag  = 'ratio';
priors_parameters_ratio.help = {['Within:between variance ratio. This ratio ' ...
    'controls the expected between-subjects ' ...
    'variability for each connection. The default is 16, meaning we expect ' ...
    'the variability in connection strengths across subjects to be 1/16 of ' ...
    'our uncertainty about connection strengths at the first level.'] '' ...
    ['Specifically, the diagonal of each precision component matrix (M.pC) ' ...
    'is set to the prior variance of the second-level parameters (M.bC) ' ...
    'divided by the ratio set here (M.beta).']};
priors_parameters_ratio.strtype = 'r';
priors_parameters_ratio.num     = [1 1];
priors_parameters_ratio.val     = {16};


% ---------------------------------------------------------------------
% priors_parameters Priors on log precision branch
% ---------------------------------------------------------------------
priors_between         = cfg_branch;
priors_between.tag     = 'priors_between';
priors_between.name    = 'Between-subjects variability';
priors_between.val     = { components ...
                           priors_parameters_ratio ...
                           priors_log_precision_mu ...
                                 priors_log_precision_var};
priors_between.help    = {['Between-subjects variability over second-' ...
     'level parameters.'], '' ...
     ['A multi-component model is used. Each component is a ' ...
      '[p x p] precision matrix given p DCM parameters, where elements on ' ...
      'the diagonal represent the precision (inverse variance) across ' ...
      'subjects of each DCM connection. These precisions are set via the ' ...
      'Within:between ratio, below. Each precision component is scaled by ' ...
      'a hyper-parameter, which is estimated from the data. The prior ' ...
      'expectation and uncertainty of these hyper-parameters are also '...
      'set below.']};

% ---------------------------------------------------------------------
% show_review Select whether to review results
% ---------------------------------------------------------------------
show_review  = cfg_menu;
show_review.tag    = 'show_review';
show_review.name   = 'Review PEB results';
show_review.labels = {'Yes','No'};
show_review.values = {1,0};
show_review.val    = {1};

% =========================================================================
% PEB specification batch
% =========================================================================

% -------------------------------------------------------------------------
% specify Specify PEB model
% -------------------------------------------------------------------------
specify      = cfg_exbranch;
specify.tag  = 'peb_specify';
specify.name = 'Specify / Estimate PEB';
specify.val  = { name model_space_mat covariates fields ...
                 bmc priors_between show_review };
specify.help = {'Specifies and estimates a second-level DCM (PEB) model.'};
specify.prog = @spm_run_create_peb;

% =========================================================================
% PEB review batch
% =========================================================================

% -------------------------------------------------------------------------
% review Review PEB model
% -------------------------------------------------------------------------
model_space_mat_op = model_space_mat;
model_space_mat_op.num = [0 Inf];
model_space_mat_op.val = {''};

review      = cfg_exbranch;
review.tag  = 'peb_review';
review.name = 'Review PEB';
review.val  = { pebmat model_space_mat_op };
review.help = {'Reviews PEB results'};
review.prog = @spm_run_dcm_peb_review;

% =========================================================================
% second_level Second level DCM batch
% =========================================================================
second_level         = cfg_choice; 
second_level.tag     = 'PEB';
second_level.name    = 'Second level';
second_level.help    = {'Parametric Empirical Bayes for DCM'};
second_level.values  = { specify review };

%==========================================================================
function out = spm_run_dcm_peb_review(job)
%==========================================================================
P   = job.pebmat;
DCM = job.model_space_mat;
spm_dcm_peb_review(P,DCM);
out = job.pebmat;

%==========================================================================
function out = spm_run_create_peb(job)
%==========================================================================

% Load provided group DCM
model_space = load(job.model_space_mat{1});
if ~isfield(model_space,'GCM')
    error('Provided file is not a valid model space');
end
GCM = model_space.GCM;

ns = size(GCM,1);

if ~isfield(GCM{1},'Ep')
    error('Please estimate DCMs before second-level analysis');
end

% DCM field(s)
if isfield(job.fields,'field_default')
    field = {'A','B'};
elseif isfield(job.fields,'field_all')
    field = 'all';
else
    field = job.fields.field_entry;
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
        error('Please ensure design matrix has one row per subject');
    end
    
    X = [X x];
    
    Xnames = [Xnames job.cov.design_mtx.cov_name];
elseif isfield(job.cov, 'regressor')
    % Design matrix entered per-regressor
    
    regressors = job.cov.regressor;
       
    for r = 1:length(regressors)
        regressor = regressors(r).cov_val;
        name      = regressors(r).cov_name;
        
        if size(regressor,1) ~= ns
            error('Please ensure regressor %d has one row per subject',r);
        end
        
        X = [X regressor];
        Xnames = [Xnames name];
    end
end

% Ensure a mean column wasn't entered accidently
if size(X,2) > 1
    bad = find(~any(diff(X(:,2:end))));
    if ~isempty(bad)
        error('Please check regressor: %d\n', bad);
    end
end

if size(X,2) ~= length(Xnames)
    error('Please ensure there is one covariate name per covariate');
end
    
% Priors / covariance components
M = struct();
M.beta   = job.priors_between.ratio;
M.hE     = job.priors_between.expectation;
M.hC     = job.priors_between.var;
M.Q      = job.priors_between.components;
M.X      = X;
M.Xnames = Xnames;

% Specify / estimate
run_bmc = (job.bmc == 1);
dir_out = fileparts(job.model_space_mat{1});
name    = job.name;
    
if run_bmc
    [BMC,PEB] = spm_dcm_bmc_peb(GCM,M,field);
    
    % Write BMC
    bmc_filename = fullfile(dir_out,['BMC_' name '.mat']);
    save(bmc_filename,'BMC');    
else
    PEB = spm_dcm_peb(GCM,M,field);
end

% Write PEB
peb_filename = fullfile(dir_out,['PEB_' name '.mat']);
save(peb_filename,'PEB');

% Review PEB
if job.show_review == 1
    spm_dcm_peb_review(peb_filename,GCM);
end
out = peb_filename;