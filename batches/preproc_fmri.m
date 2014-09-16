% preproc_fmri - preprocessing for functional MRI
%
% Need to specify (interactively, before the batch is created)
% * Whether (and if so, when) to perform slice time correction (STC)
% * The number of functional imaging sessions (note that adding/deleting a
%   session once the batch is created will break the session dependencies)
%
% Need to specify (in the created batch)
% * The slice timing acquisition details (if you're unsure, don't use STC!)
% * The functional images for each session (only in the first module)
% * The structural image (in the Segment module)
%
% PLEASE NOTE: 
% The created batch assumes structural images for different subjects will
% be in different directories (and hence that it can safely create a
% bias-corrected and skull-stripped version called simply "brain". If you
% have multiple structurals in one directory, you will need to replace
% "brain" with a unique filename for each subject.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: preproc_fmri.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('preproc_fmri', 'file') && ~isdeployed
    help preproc_fmri
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
matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
ngaus  = [1 1 2 3 4 2];
native = [1 1 1 0 0 0];
for c = 1:6 % tissue class c
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
        fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
end
matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];

%% Get directory of bias corrected image to use for brain image
iDir = length(matlabbatch) + 1;
matlabbatch{iDir}.cfg_basicio.file_dir.cfg_fileparts.files(1) = ...
    cfg_dep('Segment: Bias Corrected (1)', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));

%% Image Calculator - Create brain image (skull-stripped bias corrected)
iIC = length(matlabbatch) + 1;
for c = 1:3 % tissue class c
    matlabbatch{iIC}.spm.util.imcalc.input(c) = ...
        cfg_dep(sprintf('Segment: c%d Images', c), ...
        substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','tiss', '()',{c}, '.','c', '()',{':'}));
end
matlabbatch{iIC}.spm.util.imcalc.input(4) = ...
    cfg_dep('Segment: Bias Corrected (1)', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{iIC}.spm.util.imcalc.output = 'Brain';
matlabbatch{iIC}.spm.util.imcalc.outdir(1) = ...
    cfg_dep('Get Pathnames: Directories (unique)', ...
    substruct('.','val', '{}',{iDir}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','up'));
matlabbatch{iIC}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';

%% Coregister mean functional to Brain
iCR = length(matlabbatch) + 1;
matlabbatch{iCR}.spm.spatial.coreg.estimate.ref(1) = ...
    cfg_dep('Image Calculator: Imcalc Computed Image', ...
    substruct('.','val', '{}',{iIC}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
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

%% Normalise Brain image
iNB = length(matlabbatch) + 1;
matlabbatch{iNB}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','fordef', '()',{':'}));
matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample(1) = ...
    cfg_dep('Image Calculator: Imcalc Computed Image', ...
    substruct('.','val', '{}',{iIC}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{iNB}.spm.spatial.normalise.write.woptions.vox = [1 1 1];

%% If run as script, open matlabbatch GUI
% NB mfilename is cfg_load_jobs if loaded from matlabbatch GUI
if strcmp(mfilename, 'preproc_fmri')
    spm_jobman('interactive', matlabbatch)
end
