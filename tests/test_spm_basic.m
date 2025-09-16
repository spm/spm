function tests = test_spm_basic
% Basic unit tests for core SPM functionality (standalone compatible)
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging

% This function is designed to work with both the full MATLAB Unit Testing
% Framework and the standalone test runner

if exist('functiontests', 'builtin')
    % Full MATLAB environment
    tests = functiontests(localfunctions);
else
    % Standalone environment - return empty, tests will be discovered by parsing
    tests = [];
end

end

%==========================================================================
function test_spm_dir(testCase)
% Test that spm('Dir') returns a valid directory

d = spm('Dir');
testCase.verifyThat(d, @ischar);
testCase.verifyTrue(exist(d, 'dir') == 7);

end

%==========================================================================
function test_spm_version(testCase)
% Test that spm('Ver') returns version information

v = spm('Ver');
testCase.verifyThat(v, @ischar);
testCase.verifyTrue(~isempty(v));

end

%==========================================================================
function test_spm_memory(testCase)
% Test that spm('Memory') returns numeric values

mem = spm('Memory');
testCase.verifyThat(mem, @isnumeric);
testCase.verifyTrue(mem >= 0);

end

%==========================================================================
function test_spm_time(testCase)
% Test that spm('Time') returns a time string

t = spm('Time');
testCase.verifyThat(t, @ischar);
testCase.verifyTrue(~isempty(t));

end

%==========================================================================
function test_spm_defaults(testCase)
% Test that SPM defaults can be set and retrieved

spm('defaults', 'fmri');
% If we get here without error, the test passes
testCase.verifyTrue(true);

end
