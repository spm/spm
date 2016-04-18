function tests = test_spm_cfg_dcm_peb
% Unit Tests for spm_cfg_dcm_peb (PEB batch)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_cfg_dcm_peb.m 6770 2016-04-18 09:57:44Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function setup(testCase)
% Delete artefacts before each test

model_dir = fullfile(get_data_path(),'models');
expected_output = {fullfile(model_dir,'PEB_test.mat');
                   fullfile(model_dir,'BMA_PEB_test.mat');
                   fullfile(model_dir,'LOO_test.mat');
                   fullfile(model_dir,'GCM_simulated_abridged.mat')};

for i = 1:length(expected_output)
    if exist(expected_output{i},'file')
        spm_unlink(expected_output{i});
    end
end

% -------------------------------------------------------------------------
function test_specify_peb(testCase)
% Test specifying & estimating a PEB model

% Prepare paths
data_path = get_data_path();
model_dir = fullfile(data_path,'models');
X_file    = fullfile(data_path,'design_matrix.mat');
GCM_file  = fullfile(model_dir,'GCM_simulated.mat');

% Load design matrix
X = load(X_file);
X = X.X;

% PEB settings (different to default for test purposes)
hE = 0.1;
hC = 0.07;

% Run
clear matlabbatch;
matlabbatch{1}.spm.dcm.peb.specify.name = 'test';
matlabbatch{1}.spm.dcm.peb.specify.model_space_mat = cellstr(GCM_file);
matlabbatch{1}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).name = 'c1';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).value = X(:,1);
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).name = 'c2';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).value = X(:,2);
matlabbatch{1}.spm.dcm.peb.specify.fields.default = {'A' 'B'}';
matlabbatch{1}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.expectation = hE;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.var = hC;
matlabbatch{1}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{1}.spm.dcm.peb.specify.show_review = 0;
spm_jobman('run',matlabbatch);

% Check PEB created
expected_output = fullfile(model_dir,'PEB_test.mat');
testCase.assertTrue(exist(expected_output,'file') > 0);

% Check PEB received the right data
PEB = load(expected_output);
PEB = PEB.PEB;

expected_covariates = 3;
expected_subjects   = 30;

testCase.assertEqual(size(PEB.Ep,2), expected_covariates);
testCase.assertEqual(size(PEB.Snames,1), expected_subjects);
testCase.assertEqual(PEB.M.hE, hE);
testCase.assertEqual(PEB.M.hC, hC);

% -------------------------------------------------------------------------
function test_compare_models(testCase)

% Prepare paths
data_path = get_data_path();
model_dir = fullfile(data_path,'models');
X_file    = fullfile(data_path,'design_matrix.mat');
GCM_file  = fullfile(model_dir,'GCM_simulated.mat');

% Load design matrix
X = load(X_file);
X = X.X;

% PEB settings
hE = 0;
hC = 1/16;

% PEB batch
clear matlabbatch;
matlabbatch{1}.spm.dcm.peb.specify.name = 'test';
matlabbatch{1}.spm.dcm.peb.specify.model_space_mat = cellstr(GCM_file);
matlabbatch{1}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).name = 'c1';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).value = X(:,1);
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).name = 'c2';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).value = X(:,2);
matlabbatch{1}.spm.dcm.peb.specify.fields.default = {'A' 'B'}';
matlabbatch{1}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.expectation = hE;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.var = hC;
matlabbatch{1}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{1}.spm.dcm.peb.specify.show_review = 0;

% Model comparison batch
matlabbatch{2}.spm.dcm.peb.compare.peb_mat(1) = ...
    cfg_dep('Specify / Estimate PEB: PEB mat File(s)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', ...
    '{}',{1}, '.','val', '{}',{1}), substruct('.','peb_mat'));
matlabbatch{2}.spm.dcm.peb.compare.model_space_mat = cellstr(GCM_file);
matlabbatch{2}.spm.dcm.peb.compare.show_review = 0;

% Run
spm_jobman('run',matlabbatch);

% Check BMA created
expected_output = fullfile(model_dir,'BMA_PEB_test.mat');
testCase.assertTrue(exist(expected_output,'file') > 0);

% -------------------------------------------------------------------------
function test_search_reduced_models(testCase)

% Prepare paths
data_path = get_data_path();
model_dir = fullfile(data_path,'models');
X_file    = fullfile(data_path,'design_matrix.mat');
GCM_file  = fullfile(model_dir,'GCM_simulated.mat');

% Load design matrix
X = load(X_file);
X = X.X;

% PEB settings
hE = 0;
hC = 1/16;

% PEB batch
clear matlabbatch;
matlabbatch{1}.spm.dcm.peb.specify.name = 'test';
matlabbatch{1}.spm.dcm.peb.specify.model_space_mat = cellstr(GCM_file);
matlabbatch{1}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).name = 'c1';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(1).value = X(:,1);
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).name = 'c2';
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor(2).value = X(:,2);
matlabbatch{1}.spm.dcm.peb.specify.fields.default = {'A' 'B'}';
matlabbatch{1}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.expectation = hE;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.var = hC;
matlabbatch{1}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{1}.spm.dcm.peb.specify.show_review = 0;

% Model comparison batch
matlabbatch{2}.spm.dcm.peb.reduce_all.peb_mat(1) = ...
    cfg_dep('Specify / Estimate PEB: PEB mat File(s)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','peb_mat'));
matlabbatch{2}.spm.dcm.peb.reduce_all.model_space_mat = cellstr(GCM_file);
matlabbatch{2}.spm.dcm.peb.reduce_all.show_review = 0;

% Run
spm_jobman('run',matlabbatch);

% Check BMA created
expected_output = fullfile(model_dir,'BMA_PEB_test.mat');
testCase.assertTrue(exist(expected_output,'file') > 0);

% -------------------------------------------------------------------------
function test_loo(testCase)

% Prepare paths
data_path = get_data_path();
model_dir = fullfile(data_path,'models');
X_file    = fullfile(data_path,'design_matrix.mat');
GCM_file  = fullfile(model_dir,'GCM_simulated.mat');

% Abridge GCM
s = 11:18;
GCM_file_abridged = fullfile(model_dir,'GCM_simulated_abridged.mat');
GCM = load(GCM_file);
GCM = GCM.GCM;
GCM = GCM(s,:);
save(GCM_file_abridged,'GCM');

% Load design matrix
X = load(X_file);
X = X.X;

% PEB settings
hE = 0;
hC = 1/16;

clear matlabbatch;
matlabbatch{1}.spm.dcm.peb.predict.name = 'test';
matlabbatch{1}.spm.dcm.peb.predict.model_space_mat = cellstr(GCM_file_abridged);
matlabbatch{1}.spm.dcm.peb.predict.dcm.index = 2;
matlabbatch{1}.spm.dcm.peb.predict.cov.regressor(1).name = 'c1';
matlabbatch{1}.spm.dcm.peb.predict.cov.regressor(1).value = X(s,1);
matlabbatch{1}.spm.dcm.peb.predict.fields.custom = {'B(2,1,2)'};
matlabbatch{1}.spm.dcm.peb.predict.priors_between.ratio = 16;
matlabbatch{1}.spm.dcm.peb.predict.priors_between.expectation = hE;
matlabbatch{1}.spm.dcm.peb.predict.priors_between.var = hC;
matlabbatch{1}.spm.dcm.peb.predict.priors_glm.group_ratio = 1;

% Run
spm_jobman('run',matlabbatch);

% Check LOO created
expected_output = fullfile(model_dir,'LOO_test.mat');
testCase.assertTrue(exist(expected_output,'file') > 0);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');