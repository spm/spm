% preproc_fmri_simplified - preprocessing for functional MRI
%
% Need to specify:
% * The number of sessions (interactively, before batch is created)
% * The functional images for each session (in the "Realign & Unwarp" 
%   module of the created batch)
% * The structural image volume (in the "Segment" module)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: preproc_fmri_simplified.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('preproc_fmri_simplified', 'file') && ~isdeployed
    help preproc_fmri_simplified
end
matlabbatch = {};

%% Input: number of sessions
nSess = spm_input('Number of sessions', 1, 'n', '1', 1);

%% Realignment
iRU = 1;
for s = 1:nSess
    matlabbatch{iRU}.spm.spatial.realignunwarp.data(s).scans = ...
        '<UNDEFINED>'; %#ok<*SAGROW>
    sessdep(s) = cfg_dep(sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', s), ...
        substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','sess', '()',{s}, '.','uwrfiles'));
end

%% Segmentation
iSeg = length(matlabbatch) + 1;
matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
ngaus  = [1 1 2 3 4 2];
for c = 1:6 % tissue class c
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
        fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [0 0];
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
end
matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];

%% Coregister mean functional to bias-corrected structural
iCR = length(matlabbatch) + 1;
matlabbatch{iCR}.spm.spatial.coreg.estimate.ref(1) = ...
    cfg_dep('Segment: Bias Corrected (1)', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{iCR}.spm.spatial.coreg.estimate.source(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
    substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','meanuwr'));
    matlabbatch{iCR}.spm.spatial.coreg.estimate.other = sessdep;

%% Normalise functional images
iNF = length(matlabbatch) + 1;
matlabbatch{iNF}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','fordef', '()',{':'}));
matlabbatch{iNF}.spm.spatial.normalise.write.subj.resample = sessdep;

%% Smooth normalised functionals
iSm = length(matlabbatch) + 1;
matlabbatch{iSm}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}',{iNF}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('()',{1}, '.','files'));

%% If run as script, open matlabbatch GUI
% NB mfilename is cfg_load_jobs if loaded from matlabbatch GUI
if strcmp(mfilename, 'preproc_fmri_simplified')
    spm_jobman('interactive', matlabbatch)
end
