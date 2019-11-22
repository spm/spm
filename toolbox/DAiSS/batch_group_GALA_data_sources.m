%-----------------------------------------------------------------------
% Job saved on 11-Dec-2015 13:01:41 by cfg_util (rev $Rev: 7703 $)
% spm SPM - SPM12 (12.1)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.beamforming.data.dir = '<UNDEFINED>';
matlabbatch{1}.spm.tools.beamforming.data.D = '<UNDEFINED>';
matlabbatch{1}.spm.tools.beamforming.data.val = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.orient = 'original';
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.fdownsample = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.symmetric = 'no';
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.flip = false;
matlabbatch{2}.spm.tools.beamforming.sources.visualise = 1;
