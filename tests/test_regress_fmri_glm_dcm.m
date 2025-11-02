function tests = test_regress_fmri_glm_dcm
% Regression test for GLM and DCM for fMRI including timseries extraction
% % This script analyses the Attention to Visual Motion fMRI dataset
% available from the SPM website using DCM:
%   http://www.fil.ion.ucl.ac.uk/spm/data/attention/
% as described in the SPM docs website:
%   https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/dcm/dcm_fmri_first_level_gui/
%__________________________________________________________________________

% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

tests = functiontests({@setup;
                       @test_glm_parametric_volterra;
                       @test_regress_glm_dcm});

function setup(testCase)

% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

% -------------------------------------------------------------------------
function test_glm_parametric_volterra(testCase)
% Test GLM specification with extra features: 
% parametric regressors, volterra kernels and 2 basis functions

data_path = fullfile(spm('Dir'),'tests','data','attention');

options.bf = 2;
options.volterra = true;
options.pmod = true;

% Specify and estimate GLM
spm_setup_data_attention(data_path,options);

% We expect 30 design matrix columns (9 first order, 21 second order)
SPM = load(fullfile(data_path,'GLM','SPM.mat'));
SPM = SPM.SPM;
testCase.verifyTrue(size(SPM.xX.X,2)==30);

% -------------------------------------------------------------------------
function test_regress_glm_dcm(testCase)
% Test a simpler GLM with timeseries extraction and DCM

data_path = fullfile(spm('Dir'),'tests','data','attention');

% Specify and estimate GLM, extract timeseries and run DCM
options.dcm = true;
spm_setup_data_attention(data_path,options);

% Load estimated DCM and BMR results
results = load(fullfile(data_path,'GCM_test.mat'));
GCM  = results.GCM;
post = results.post;

% Get the explained variance
DCM = spm_dcm_fmri_check(GCM{1},true);
exp_var = DCM.diagnostics(1);
max_connection = DCM.diagnostics(2);

% Check the forward model was the winner (Pp > 95%)
testCase.verifyTrue(post(1) > 0.95); 

% Check the explained variance was high (> 75%)
testCase.verifyTrue(exp_var > 75);

% Check the largest neural connection was non-trivial (> 0.5Hz)
testCase.verifyTrue(max_connection > 0.5);