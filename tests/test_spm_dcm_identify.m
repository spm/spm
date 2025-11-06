classdef test_spm_dcm_identify < matlab.unittest.TestCase
% Unit Tests for test_spm_dcm_identify
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


methods (Test)

% -------------------------------------------------------------------------
function test_identify_dcm_fmri(testCase)

data_path = test_spm_dcm_identify.get_data_path();

DCM = load(fullfile(data_path,'DCM_fMRI.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'fMRI');
end

% -------------------------------------------------------------------------
function test_identify_dcm_fmri_csd(testCase)

data_path = test_spm_dcm_identify.get_data_path();

DCM = load(fullfile(data_path,'DCM_fMRI_CSD.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'fMRI_CSD');
end

% -------------------------------------------------------------------------
function test_identify_erp(testCase)

data_path = test_spm_dcm_identify.get_data_path();

DCM = load(fullfile(data_path,'DCM_IMG_ERP_CMC.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'ERP');
end

% -------------------------------------------------------------------------
function test_identify_csd(testCase)

data_path = test_spm_dcm_identify.get_data_path();

DCM = load(fullfile(data_path,'DCM_LFP_CSD_CMC.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'CSD');
end

% -------------------------------------------------------------------------
function test_identify_ind(testCase)

data_path = test_spm_dcm_identify.get_data_path();

DCM = load(fullfile(data_path,'DCM_IND.mat'));
DCM = DCM.DCM;

model = spm_dcm_identify(DCM);

testCase.assertEqual(model,'IND');
end

end % methods (Test)

methods (Static, Access = private)
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_identify');
end
end % methods (Static, Access = private)

end % classdef