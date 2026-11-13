%% Test SPM for fMRI

% Path to test data
data_path = fullfile(spm('Dir'),'tests','data','attention');

% Clear old output if present
output = fullfile(data_path,'GCM_test.mat');
if exist(output,'file')
    spm_unlink(output);
end

% Run GLM, timseries extraction and DCM for fMRI
options = struct();
options.dcm = true;
spm_setup_data_attention(data_path,options);

% Confirm everything ran to completion
assert(exist(output,'file') > 0);