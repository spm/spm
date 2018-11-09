function tests = test_spm_dcm_fmri_csd

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_endogenous(testCase)
% Tests a model with no driving inputs

spm('defaults','fmri');
DCM = load(fullfile(get_data_path(), 'DCM_attention_CSD_endogenous_r7259.mat'));
DCM = DCM.DCM;

DCM = spm_dcm_fmri_csd(DCM);

testCase.assertTrue(isfield(DCM,'Ep'));
testCase.assertTrue(isfield(DCM,'Cp'));
testCase.assertTrue(isfield(DCM,'F'));

% -------------------------------------------------------------------------
function test_driving(testCase)
% Tests a model with no driving inputs

spm('defaults','fmri');
DCM = load(fullfile(get_data_path(), 'DCM_attention_CSD_r7259.mat'));
DCM = DCM.DCM;

DCM = spm_dcm_fmri_csd(DCM);

testCase.assertTrue(isfield(DCM,'Ep'));
testCase.assertTrue(isfield(DCM,'Cp'));
testCase.assertTrue(isfield(DCM,'F'));
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_fmri_csd');