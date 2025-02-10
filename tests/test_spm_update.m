function tests = test_spm_update
% Unit Tests for spm_update
%__________________________________________________________________________

% Copyright (C) 2017-2025 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_version_check(testCase)
% Checks whether the server can be contacted. Makes max 3 attempts.

attempts = 3;
failed = false(attempts,1);
for i = 1:attempts
    sts = spm_update();
    failed(i) = all(isnan(sts)) | all(isinf(sts)) | isempty(sts);
    
    if failed(i)
        % Failure - pause a second
        pause(1);
    else
        % Success - stop here
        break
    end
end

testCase.assertFalse(all(failed));