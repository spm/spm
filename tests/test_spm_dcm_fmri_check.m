classdef test_spm_dcm_fmri_check < matlab.unittest.TestCase
% Unit Tests for spm_dcm_fmri_check
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


methods (Test)

% -------------------------------------------------------------------------
function test_DCM(testCase)
% Simply tests the function doesn't crash with one DCM

import matlab.unittest.constraints.*

data_path = test_spm_dcm_fmri_check.get_data_path();
dcm_file  = fullfile(data_path,'DCM_s1_m1.mat');

DCM = spm_dcm_fmri_check(dcm_file, true);

testCase.fatalAssertThat(DCM, HasField('diagnostics'));
testCase.verifyThat(DCM.diagnostics, HasElementCount(3));
end

% -------------------------------------------------------------------------
function test_GCM(testCase)
% Simply tests the function doesn't crash with a GCM array

import matlab.unittest.constraints.*

data_path = test_spm_dcm_fmri_check.get_data_path();
gcm_file  = fullfile(data_path,'GCM_simulated.mat');

GCM = spm_dcm_fmri_check(gcm_file, true);

testCase.verifyThat(GCM, IsOfClass('cell'));
testCase.fatalAssertThat(GCM{1}, HasField('diagnostics'));
testCase.verifyThat(GCM{1}.diagnostics, HasElementCount(3));
end

end % methods (Test)

methods (Static, Access = private)
% -------------------------------------------------------------------------
function data_path = get_data_path()
data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region', 'models');
end
end % methods (Static, Access = private)

end % classdef