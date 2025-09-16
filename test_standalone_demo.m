% Simple test script to demonstrate standalone compilation approach
% This would be used to create a minimal standalone test

fprintf('=== SPM Standalone Test Demo ===\n');

% Add SPM to path
addpath(fileparts(mfilename('fullpath')));

% Initialize SPM
try
    spm('defaults', 'fmri');
    spm_get_defaults('cmdline', true);
    fprintf('✓ SPM initialized successfully\n');
catch ME
    fprintf('✗ SPM initialization failed: %s\n', ME.message);
    return;
end

% Test basic SPM functionality
try
    d = spm('Dir');
    fprintf('✓ SPM directory: %s\n', d);
    
    v = spm('Ver');
    fprintf('✓ SPM version: %s\n', v);
    
    if exist('runtests', 'builtin')
        fprintf('✓ runtests function is available\n');
        
        % Try to run a simple test
        test_file = fullfile(spm('Dir'), 'tests', 'test_spm.m');
        if exist(test_file, 'file')
            fprintf('  Attempting to run test_spm.m...\n');
            results = runtests(test_file);
            
            if ~isempty(results)
                fprintf('  ✓ Test completed: %d passed, %d failed\n', ...
                    sum([results.Passed]), sum([results.Failed]));
            else
                fprintf('  ✗ No test results returned\n');
            end
        else
            fprintf('  - test_spm.m not found\n');
        end
    else
        fprintf('✗ runtests function not available in this environment\n');
    end
    
catch ME
    fprintf('✗ Error during testing: %s\n', ME.message);
end

fprintf('=== Test completed ===\n');
