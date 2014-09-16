% stats2_between_subject - second-level fMRI or PET/VBM statistics
%
% Need to specify:
% * Factorial design specification: Directory
% * Factorial design specification: Scans
% * Contrast Manager: one or more contrasts
%
% Consider specifying:
% * The details of your design (!)
% * An explicit mask (for VBM)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: stats2_between_subject.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('stats2_between_subject', 'file') && ~isdeployed
    help stats2_between_subject
end
matlabbatch = {};

%%
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = '<UNDEFINED>';

matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
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
if strcmp(mfilename, 'stats2_between_subject')
    spm_jobman('interactive', matlabbatch)
end
