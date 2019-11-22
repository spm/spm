%-----------------------------------------------------------------------
% Job saved on 04-Sep-2015 11:55:24 by cfg_util (rev $Rev: 7703 $)
% spm SPM - SPM12 (12.1)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------#
matlabbatch{1}.spm.tools.beamforming.data.val = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{1}.spm.tools.beamforming.data.overwrite = 0;
matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.orient = 'unoriented';
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.fdownsample = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.symmetric = 'no';
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.flip = false;
matlabbatch{2}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{3}.spm.tools.beamforming.features.woi = [-100 0
                                                     100 200];
matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{3}.spm.tools.beamforming.features.plugin.vbfa.nl = 5;
matlabbatch{3}.spm.tools.beamforming.features.plugin.vbfa.nem = 50;
matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 0;
matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.champagne.nem = 100;
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.champagne.vcs = 0;
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.champagne.nupd = 0;
matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.woi = [100 200];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.foi = [0 Inf];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.result = 'bycondition';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.scale = 0;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.write.plugin.gifti.normalise = 'all';
matlabbatch{6}.spm.tools.beamforming.write.plugin.gifti.space = 'mni';
matlabbatch{6}.spm.tools.beamforming.write.plugin.gifti.visualise = 2;
