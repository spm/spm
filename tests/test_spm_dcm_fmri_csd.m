classdef test_spm_dcm_fmri_csd < matlab.unittest.TestCase
% Unit Tests for spm_dcm_fmri_csd
%__________________________________________________________________________

% Copyright (C) 2017-2022 Wellcome Centre for Human Neuroimaging

methods (TestClassSetup)
    function setupSPM(testCase)
        spm_get_defaults('dcm.verbose',false);
        spm_get_defaults('cmdline',true);
    end
end % methods (TestClassSetup)


methods (Test)

% -------------------------------------------------------------------------
function test_endogenous(testCase)
% Tests a model with no driving inputs
DCM = load(fullfile(test_spm_dcm_fmri_csd.get_data_path(), 'DCM_attention_CSD_endogenous_r7259.mat'));
DCM = DCM.DCM;

DCM = spm_dcm_fmri_csd(DCM);

testCase.assertTrue(isfield(DCM,'Ep'));
testCase.assertTrue(isfield(DCM,'Cp'));
testCase.assertTrue(isfield(DCM,'F'));
end

% -------------------------------------------------------------------------
function test_driving(testCase)
% Tests a model with no driving inputs

DCM = load(fullfile(test_spm_dcm_fmri_csd.get_data_path(), 'DCM_attention_CSD_r7259.mat'));
DCM = DCM.DCM;

DCM = spm_dcm_fmri_csd(DCM);

testCase.assertTrue(isfield(DCM,'Ep'));
testCase.assertTrue(isfield(DCM,'Cp'));
testCase.assertTrue(isfield(DCM,'F'));
end

end % methods (Test)

methods (Static, Access = private)
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_fmri_csd');
end
end % methods (Static, Access = private)

end % classdef