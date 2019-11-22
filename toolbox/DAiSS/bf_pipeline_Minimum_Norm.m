%-----------------------------------------------------------------------
% Job saved on 02-Apr-2014 13:13:09 by cfg_util (rev $Rev: 7703 $)
% spm SPM - SPM12b (beta)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.source.headmodel.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
matlabbatch{2}.spm.tools.beamforming.data.dir = '<UNDEFINED>';
matlabbatch{2}.spm.tools.beamforming.data.D(1) = cfg_dep('M/EEG head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{2}.spm.tools.beamforming.data.val = 1;
matlabbatch{2}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{2}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{2}.spm.tools.beamforming.data.overwrite = 1;
matlabbatch{3}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{3}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.mesh.orient = 'original';
matlabbatch{3}.spm.tools.beamforming.sources.plugin.mesh.fdownsample = 1;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.mesh.flip = false;
matlabbatch{3}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{4}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{4}.spm.tools.beamforming.features.woi = [-Inf Inf];
matlabbatch{4}.spm.tools.beamforming.features.modality = {'MEG'};
matlabbatch{4}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.foi = [0 Inf];
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
matlabbatch{4}.spm.tools.beamforming.features.regularisation.manual.lambda = 0;
matlabbatch{4}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{5}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.minimumnorm.snr = 5;
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.minimumnorm.trunc = 0;
matlabbatch{6}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.woi = [-Inf Inf];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.scale = 0;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
matlabbatch{7}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{7}.spm.tools.beamforming.write.plugin.gifti.normalise = 'no';
matlabbatch{7}.spm.tools.beamforming.write.plugin.gifti.space = 'mni';
matlabbatch{7}.spm.tools.beamforming.write.plugin.gifti.visualise = 2;
