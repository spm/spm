%-----------------------------------------------------------------------
% Job saved on 27-Nov-2015 18:36:22 by cfg_util (rev $Rev: 7703 $)
% spm SPM - SPM12 (12.1)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.beamforming.group.BF = '<UNDEFINED>';
matlabbatch{1}.spm.tools.beamforming.group.prefix = '';
matlabbatch{1}.spm.tools.beamforming.group.plugin.batch.batchfile = {'C:\spm12\toolbox\DAiSS\batch_group_GALA_data_sources.m'};
matlabbatch{2}.spm.tools.beamforming.group.BF(1) = cfg_dep('Group analysis: BF.mat files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.group.prefix = '';
matlabbatch{2}.spm.tools.beamforming.group.plugin.GALA.iter = 3;
matlabbatch{3}.spm.tools.beamforming.group.BF(1) = cfg_dep('Group analysis: BF.mat files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.group.prefix = '';
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.measure = 'lJcov';
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.spread = 4;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.threshold = 0.01;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.mincorr = 0.7;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.maxsize = 50;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.distratio1 = 4;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.distratio2 = 2;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.cluster.maxclust.maxclustsize = 30;
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.linkmeth = 'complete';
matlabbatch{3}.spm.tools.beamforming.group.plugin.functionalROI.similarity = 0;
matlabbatch{4}.spm.tools.beamforming.group.BF(1) = cfg_dep('Group analysis: BF.mat files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.group.prefix = '';
matlabbatch{4}.spm.tools.beamforming.group.plugin.batch.batchfile = {'C:\spm12\toolbox\DAiSS\batch_group_GALA_write.m'};