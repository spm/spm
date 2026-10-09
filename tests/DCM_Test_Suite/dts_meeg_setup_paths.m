function root = dts_meeg_setup_paths()
% Add SPM regression suite helper/test paths needed by scripts and tests.
% Authored by Pranay Yadav in 2026

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fileparts(root));
addpath(fullfile(spm('Dir'), 'matlabbatch'));
spm_get_defaults('cmdline', true);
spm_get_defaults('dcm.verbose', false);
spm('CmdLine', true);
spm('defaults', 'eeg');
end
