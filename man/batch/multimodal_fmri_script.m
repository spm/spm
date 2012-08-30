% To use this example, please 
% 1) place this file in your MATLAB path
% 2) place the 'multimodal_fmri_template.m' in the 'fMRI' folder of your
%    copy of the multimodal dataset
% 3) edit the line PLEASE_ENTER_YOUR_PATH_TO_MULTIMODAL_DATASET below

% Inputs per subject
% * subject dir
% * EPI images
%   - session 1
%   - session 2
% * Anatomy image
% * Multiple conditions
%   - session 1
%   - session 2
% Assumed data organisation
% multimodal/fMRI/{Session1,Session2}

studydir  = PLEASE_ENTER_YOUR_PATH_TO_MULTIMODAL_DATASET;
subjdirs  = {'.'};
jobs      = cell(size(subjdirs));
[jobs{:}] = deal(fullfile(studydir,'fMRI','multimodal_fmri_template.m'));

inputs    = cell(6,numel(subjdirs));
for sub = 1:numel(subjdirs)
    inputs{1,sub} = {fullfile(studydir,subjdirs{sub})};
    inputs{2,sub} = cellstr(spm_select('FPList',fullfile(inputs{1,sub}{1}, ...
                                                  'fMRI','Session1'),'^f.*\.img'));
    inputs{3,sub} = cellstr(spm_select('FPList',fullfile(inputs{1,sub}{1}, ...
                                                  'fMRI','Session2'),'^f.*\.img'));
    inputs{4,sub} = {fullfile(inputs{1,sub}{1},'sMRI','smri.img')};
    inputs{5,sub} = {fullfile(inputs{1,sub}{1},'fMRI','trials_ses1.mat')};
    inputs{6,sub} = {fullfile(inputs{1,sub}{1},'fMRI','trials_ses2.mat')};
end
spm('defaults','fmri');
spm_jobman('run',jobs,inputs{:});
