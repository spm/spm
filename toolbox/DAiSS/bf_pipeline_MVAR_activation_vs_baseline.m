%-----------------------------------------------------------------------
% Job saved on 08-Feb-2017 12:35:41 by cfg_util (rev $Rev: 7703 $)
% spm SPM - SPM12 (12.3)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.source.headmodel.D = {''};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'FIL_CTF_L';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'FIL_CTF_R';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
matlabbatch{2}.spm.tools.beamforming.data.dir = '<UNDEFINED>';
matlabbatch{2}.spm.tools.beamforming.data.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{2}.spm.tools.beamforming.data.val = 1;
matlabbatch{2}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{2}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{2}.spm.tools.beamforming.data.overwrite = 0;
matlabbatch{3}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{3}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.grid.resolution = 10;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';
matlabbatch{3}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{4}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{4}.spm.tools.beamforming.features.woi = [-Inf Inf];
matlabbatch{4}.spm.tools.beamforming.features.modality = {'MEG'};
matlabbatch{4}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.foi = [0 100];
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
matlabbatch{4}.spm.tools.beamforming.features.regularisation.manual.lambda = 5;
matlabbatch{4}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{5}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
matlabbatch{6}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.isdesign.custom.whatconditions.all = 1;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.isdesign.custom.contrast = [-1 1];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.isdesign.custom.woi = [-100 0
                                                                                   100 200];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.datafeatures = 'sumpower';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.foi = [5 15
                                                                   15 30
                                                                   30 60
                                                                   60 90];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.result = 'BIC';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.sametrials = false;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_mv.modality = 'MEG';
matlabbatch{7}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{7}.spm.tools.beamforming.write.plugin.nifti.normalise = 'no';
matlabbatch{7}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';
