% List of open inputs
% Named Directory Selector: Directory - cfg_files
% Realign: Estimate & Reslice: Session - cfg_files
% Coreg: Estimate: Source Image - cfg_files
% fMRI model specification: Multiple conditions - cfg_files
nrun = X; % enter the number of runs here
jobfile = {fullfile(spm('dir'),'man','batch','face_single_subject_template.m')};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Named Directory Selector: Directory - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Realign: Estimate & Reslice: Session - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Coreg: Estimate: Source Image - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % fMRI model specification: Multiple conditions - cfg_files
end
spm('defaults','fmri');
spm_jobman('run',jobs,inputs{:});
