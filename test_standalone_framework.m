% Test script to verify standalone test functionality
% This simulates what would happen in the standalone version

fprintf('Testing SPM Standalone Test Framework...\n\n');

% Test 1: Basic functionality test
fprintf('=== Test 1: Basic SPM Functions ===\n');
try
    exit_code = spm_standalone_tests();
    fprintf('Basic tests exit code: %d\n\n', exit_code);
catch ME
    fprintf('Error running basic tests: %s\n\n', ME.message);
end

% Test 2: Specific test
fprintf('=== Test 2: Specific Test (spm_basic) ===\n');
try
    exit_code = spm_standalone_tests('spm_basic');
    fprintf('Specific test exit code: %d\n\n', exit_code);
catch ME
    fprintf('Error running specific test: %s\n\n', ME.message);
end

% Test 3: Test a known working test file
fprintf('=== Test 3: Platform Test ===\n');
try
    exit_code = spm_standalone_tests('spm_platform');
    fprintf('Platform test exit code: %d\n\n', exit_code);
catch ME
    fprintf('Error running platform test: %s\n\n', ME.message);
end

fprintf('Test script completed.\n');
