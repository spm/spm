function tests = test_spm_update
% Unit Tests for spm_update
%__________________________________________________________________________

% Copyright (C) 2017-2025 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_version_check(testCase)
% Checks whether the server can be contacted. Makes max 3 attempts.

attempts = 1;
failed = false(attempts,1);
for i = 1:attempts
    [sts, msg] = spm_update();
    
    failed(i) = all(isnan(sts)) | all(isinf(sts)) | isempty(sts);
    
    if failed(i) && attempts > 1
        % Failure - pause a second
        pause(1);
    else
        % Success - stop here
        break
    end
end

if all(failed)
    warning('Could not run spm_update. This is a known issue when tests are run on Github');
    disp(sts);
    disp(msg);
end