function tests = test_spm_update
% Unit Tests for spm_update
%__________________________________________________________________________

% Copyright (C) 2017-2025 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_version_check(testCase)
% Checks whether the server can be contacted
sts = spm_update();
testCase.assertFalse(all(isnan(sts)) | all(isinf(sts)) | isempty(sts));