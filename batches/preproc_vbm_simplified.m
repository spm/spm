% preproc_vbm_simplified - preprocessing for voxel-based morphometry
%
% Need to specify:
% * Structural images for all subjects
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: preproc_vbm_simplified.m 6177 2014-09-16 10:44:55Z guillaume $

if exist('preproc_vbm_simplified', 'file') && ~isdeployed
    help preproc_vbm_simplified
end
matlabbatch = {};

%% Segmentation
iSeg = 1;
ngaus  = [1 1 2 3 4 2];
native = [
    1 0 0 0 0 0
    1 1 0 0 0 0];
for c = 1:6 % tissue class c
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
        fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = native(:, c)';
    matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
end

%% 
iDar = 2;
matlabbatch{iDar}.spm.tools.dartel.warp.images{1}(1) = ...
    cfg_dep('Segment: rc1 Images', substruct('.','val', '{}',{iSeg}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{1}, '.','rc', '()',{':'}));
matlabbatch{iDar}.spm.tools.dartel.warp.images{2}(1) = ...
    cfg_dep('Segment: rc2 Images', substruct('.','val', '{}',{iSeg}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{2}, '.','rc', '()',{':'}));

%% Normalise segmentations to MNI space
iNS = 3;
matlabbatch{iNS}.spm.tools.dartel.mni_norm.template(1) = ...
    cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', ...
    substruct('.','val', '{}',{iDar}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','template', '()',{7}));
matlabbatch{iNS}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1) = ...
    cfg_dep('Run Dartel (create Templates): Flow Fields', ...
    substruct('.','val', '{}',{iDar}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files', '()',{':'}));
matlabbatch{iNS}.spm.tools.dartel.mni_norm.data.subjs.images{1}(1) = ...
    cfg_dep('Segment: c1 Images', ...
    substruct('.','val', '{}',{iSeg}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{iNS}.spm.tools.dartel.mni_norm.preserve = 1;
matlabbatch{iNS}.spm.tools.dartel.mni_norm.fwhm = [8 8 8];

%% If run as script, open matlabbatch GUI
% NB mfilename is cfg_load_jobs if loaded from matlabbatch GUI
if strcmp(mfilename, 'preproc_vbm_simplified')
    spm_jobman('interactive', matlabbatch)
end
