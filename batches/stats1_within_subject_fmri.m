% stats1_within_subject_fmri - first-level fMRI statistics
%
% Need to specify:
% * fMRI model specification: Directory
% * fMRI model specification: Units for design
% * fMRI model specification: Interscan interval
% * fMRI model specification: Scans
% * Contrast Manager: one or more contrasts
%
% Consider specifying:
% * Conditions (not needed for resting-state studies)
% * Multiple regressors - the realignment parameters from motion correction
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: stats1_within_subject_fmri.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('stats1_within_subject_fmri', 'file') && ~isdeployed
    help stats1_within_subject_fmri
end
matlabbatch = {};

%%
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = '<UNDEFINED>';

matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
    cfg_dep('fMRI model specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));

matlabbatch{3}.spm.stats.con.spmmat = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));

matlabbatch{4}.spm.stats.results.spmmat = ...
    cfg_dep('Contrast Manager: SPM.mat File', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.contrasts = Inf;

%% If run as script, open matlabbatch GUI
% NB mfilename is cfg_load_jobs if loaded from matlabbatch GUI
if strcmp(mfilename, 'stats1_within_subject_fmri')
    spm_jobman('interactive', matlabbatch)
end
