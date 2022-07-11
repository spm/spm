function tests = test_spm_dcm_identify
% Unit Tests for test_spm_dcm_identify
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_identify_dcm_fmri(testCase)

data_path = get_data_path();

DCM = load(fullfile(data_path,'DCM_fMRI.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'fMRI');

% -------------------------------------------------------------------------
function test_identify_dcm_fmri_csd(testCase)

data_path = get_data_path();

DCM = load(fullfile(data_path,'DCM_fMRI_CSD.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'fMRI_CSD');

% -------------------------------------------------------------------------
function test_identify_erp(testCase)

data_path = get_data_path();

DCM = load(fullfile(data_path,'DCM_IMG_ERP_CMC.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'ERP');

% -------------------------------------------------------------------------
function test_identify_csd(testCase)

data_path = get_data_path();

DCM = load(fullfile(data_path,'DCM_LFP_CSD_CMC.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'CSD');

% -------------------------------------------------------------------------
function test_identify_ind(testCase)

data_path = get_data_path();

DCM = load(fullfile(data_path,'DCM_IND.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'IND');

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_identify');
