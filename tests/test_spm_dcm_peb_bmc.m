function tests = test_spm_dcm_peb_bmc
% Unit Tests for test_spm_dcm_peb_bmc
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_model_search(testCase)

data_path = get_data_path();

PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;

[BMA,BMR] = spm_dcm_peb_bmc(PEB(1));

% Get BMA parameters
[n,m] = size(PEB(1).Ep);
Ep = reshape(full(BMA.Ep),[n,m]);
Pp = reshape(full(BMA.Pp),[n,m]);

effect     = 2; % group difference
connection = 5; % r1->r2

% There should be an effect of group on the connection from R1->R2 of 0.2Hz 
testCase.assertEqual(Ep(connection,effect),0.2, 'AbsTol', 0.05);
testCase.assertTrue(Pp(connection,effect) > 0.9);

% There should be no effect of group elsewhere
Pp_others = Pp; 
Pp_others(connection,:) = [];
testCase.assertTrue(all(Pp_others(connection,effect) < 0.95));

close all;

% -------------------------------------------------------------------------
function test_specific_model_comparison(testCase)

data_path = get_data_path();

% Get PEB full model
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB(1);

% Get template models
GCM_templates = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM_templates = GCM_templates.GCM;

% Run model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),GCM_templates(1,:));

% Get BMA parameters
[n,m] = size(PEB(1).Ep);
Ep = reshape(full(BMA.Ep),[n,m]);

effect     = 2; % group difference
connection = 5; % r1->r2

% There should be common effects on the forward connection only
testCase.assertTrue(BMA.Pw(1) > 0.95);
testCase.assertTrue(BMA.Pw(2) < 0.95);

% There should be an effect of group on the connection from R1->R2 of 0.2Hz 
testCase.assertEqual(Ep(connection,effect),0.2, 'AbsTol', 0.05);
testCase.assertTrue(BMA.Px(1) > 0.95);
testCase.assertTrue(BMA.Px(2) < 0.95);

close all;

% -------------------------------------------------------------------------
function test_self_connection_onoff(testCase)
% Implemented to confirm that self-connections on A in DCM for fMRI models
% can be varied across models.
data_path = get_data_path();

% Get PEB full model
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB(1);

% Get template models
GCM_templates = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM_templates = GCM_templates.GCM;

% Create model space
GCM_templates      = GCM_templates(1,1:2);
GCM_templates{2}   = rmfield(GCM_templates{2},'M');
GCM_templates{2}.b = GCM_templates{1}.b;
GCM_templates{2}.a(2,2) = 0;
GCM_templates{2}.b(2,1,2) = 0;

% Run model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),GCM_templates);

% There should be two parameters varying across models
testCase.assertTrue(length(BMA.Kname) == 2);

close all;

% -------------------------------------------------------------------------
function test_bmc_peb_of_pebs(testCase)
% Tests model comparison with hierarchical models

data_path = get_data_path();

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Prepare group level design matrix
X = load(fullfile(data_path,'design_matrix.mat'));
X = X.X;
ns = size(X,1);
X  = [ones(ns,1) X];

% Assign subjects to groups
g1 = (X(:,2) == -1);
g2 = (X(:,2) == 1);

% Get a random covariate (e.g. age)
age_g1 = X(g1,3) - mean(X(g1,3));
age_g2 = X(g2,3) - mean(X(g2,3));

% Estimate PEB
M = struct();
M.Q = 'all';
M.Xnames = {'Mean','Age'};

M.X = [ones(ns/2,1) age_g1];
PEB1 = spm_dcm_peb(GCM(g1,1),M);

M.X = [ones(ns/2,1) age_g2];
PEB2 = spm_dcm_peb(GCM(g2,1),M);

% Group model
M        = struct();
M.X      = [1 1; 1 -1];
M.Q      = 'none';
M.Xnames = {'Commonalities','Differences'};
PEB      = spm_dcm_peb({PEB1;PEB2},M);

% Get template models
GCM_templates = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM_templates = GCM_templates.GCM;

% Create model space
GCM_templates      = GCM_templates(1,1:2);
GCM_templates{1}   = rmfield(GCM_templates{1},'M');
GCM_templates{2}   = rmfield(GCM_templates{2},'M');
GCM_templates{2}.b = GCM_templates{1}.b;
GCM_templates{2}.b(2,1,2) = 0;

% Run model comparison: we expect one parameter to vary
[BMA,BMR] = spm_dcm_peb_bmc(PEB,GCM_templates,'Mean');
assert(length(BMA.Kname) == 1);

close all

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');
