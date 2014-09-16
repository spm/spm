% preproc_fmri_without_anat - preprocessing for fMRI without anatomical MRI
%
% Need to specify (interactively, before the batch is created)
% * Whether (and if so, when) to perform slice time correction (STC)
% * The number of functional imaging sessions (note that adding/deleting a
%   session once the batch is created will break the session dependencies)
%
% Need to specify (in the created batch)
% * The slice timing acquisition details (if you're unsure, don't use STC!)
% * The functional images for each session (only in the first module)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: preproc_fmri_without_anat.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('preproc_fmri_without_anat', 'file') && ~isdeployed
    help preproc_fmri_without_anat
end
matlabbatch = {};

%% Inputs: slice time correction option; number of sessions
STC   = spm_input('Slice time correction', 1, 'm', {
    'Do not perform STC'
    'STC before realignment'
    'STC after realignment'
    }, [0 1 2], 1);
nSess = spm_input('Number of sessions', '+1', 'n', '1', 1);

%% Slice time correction and/or Realign & Unwarp
switch STC
    case 0
        iRU = 1;
        for s = 1:nSess
            matlabbatch{iRU}.spm.spatial.realignunwarp.data(s).scans = ...
                '<UNDEFINED>'; %#ok<*SAGROW>
            sessdep(s) = cfg_dep(sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', s), ...
                substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','sess', '()',{s}, '.','uwrfiles'));
        end
    case 1
        iST = 1; iRU = 2;
        for s = 1:nSess % session s
            matlabbatch{iST}.spm.temporal.st.scans{1, s} = '<UNDEFINED>';
            matlabbatch{iRU}.spm.spatial.realignunwarp.data(s).scans(1) = ...
                cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)', s), ...
                substruct('.','val', '{}',{iST}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{s}, '.','files'));
            sessdep(s) = cfg_dep(sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', s), ...
                substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','sess', '()',{s}, '.','uwrfiles'));
        end
    case 2
        iRU = 1; iST = 2;
        for s = 1:nSess
            matlabbatch{iRU}.spm.spatial.realignunwarp.data(s).scans = '<UNDEFINED>';
            matlabbatch{iST}.spm.temporal.st.scans{1, s} = cfg_dep(...
                sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', s), ...
                substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','sess', '()',{s}, '.','uwrfiles'));
            sessdep(s) = cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)', s), ...
                substruct('.','val', '{}',{iST}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{s}, '.','files'));
        end
end

%% Segmentation
iSeg = length(matlabbatch) + 1;
matlabbatch{iSeg}.spm.spatial.preproc.channel.vols(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
    substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','meanuwr'));
ngaus  = [1 1 2 2 2 2]; 
% Note: The number of Gaussians could probably be improved; it might also
% be better to use a single TPM for GM/WM, due to the poor contrast in T2*
native = [1 1 1 0 0 0]; % not needed, but helpful to check quality of seg.
for c = 1:6 % tissue class c
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
        fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
end
matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];

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

%% Normalise mean functional
iNM = length(matlabbatch) + 1;
matlabbatch{iNM}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','fordef', '()',{':'}));
matlabbatch{iNM}.spm.spatial.normalise.write.subj.resample(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
    substruct('.','val', '{}',{iRU}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','meanuwr'));
matlabbatch{iNM}.spm.spatial.normalise.write.woptions.vox = [1 1 1];

%% If run as script, open matlabbatch GUI
% NB mfilename is cfg_load_jobs if loaded from matlabbatch GUI
if strcmp(mfilename, 'preproc_fmri_without_anat')
    spm_jobman('interactive', matlabbatch)
end
