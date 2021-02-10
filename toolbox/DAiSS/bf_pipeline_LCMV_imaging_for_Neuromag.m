%-----------------------------------------------------------------------
% Job saved on 10-Feb-2021 12:40:00 by cfg_util (rev $Rev: 8061 $)
% cfg_basicio BasicIO - Unknown
% dtijobs DTI tools - Unknown
% impexp_NiftiMrStruct NiftiMrStruct - Unknown
% spm SPM - SPM12 (12.6)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.preproc.reduce.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.preproc.reduce.timewin = [-Inf Inf];
matlabbatch{1}.spm.meeg.preproc.reduce.channels{1}.type = 'MEG';
matlabbatch{1}.spm.meeg.preproc.reduce.channels{2}.type = 'MEGPLANAR';
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.Dnoise = {''};
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.method = 'tikhonov';
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.separate = true;
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.lambda = 0;
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.lambda_percent = false;
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.kappa = 60;
matlabbatch{1}.spm.meeg.preproc.reduce.method.whiten.tolerance = 0;
matlabbatch{1}.spm.meeg.preproc.reduce.keeporig = false;
matlabbatch{1}.spm.meeg.preproc.reduce.keepothers = false;
matlabbatch{1}.spm.meeg.preproc.reduce.prefix = 'R';
matlabbatch{2}.spm.meeg.source.headmodel.D(1) = cfg_dep('Data reduction: Reduced M/EEG dataset', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.source.headmodel.val = 1;
matlabbatch{2}.spm.meeg.source.headmodel.comment = '';
matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
matlabbatch{2}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{2}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
matlabbatch{3}.spm.tools.beamforming.data.dir = '<UNDEFINED>';
matlabbatch{3}.spm.tools.beamforming.data.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{3}.spm.tools.beamforming.data.val = 1;
matlabbatch{3}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{3}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{3}.spm.tools.beamforming.data.overwrite = 1;
matlabbatch{4}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{4}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{4}.spm.tools.beamforming.sources.plugin.grid.resolution = 10;
matlabbatch{4}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';
matlabbatch{4}.spm.tools.beamforming.sources.plugin.grid.constrain = 'iskull';
matlabbatch{4}.spm.tools.beamforming.sources.normalise_lf = false;
matlabbatch{4}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{5}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{5}.spm.tools.beamforming.features.woi = '<UNDEFINED>';
matlabbatch{5}.spm.tools.beamforming.features.modality = {
                                                          'MEG'
                                                          'MEGPLANAR'
                                                          }';
matlabbatch{5}.spm.tools.beamforming.features.fuse = 'meg';
matlabbatch{5}.spm.tools.beamforming.features.cross_terms = 'all';
matlabbatch{5}.spm.tools.beamforming.features.plugin.cov.foi = [0 Inf];
matlabbatch{5}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
matlabbatch{5}.spm.tools.beamforming.features.regularisation.mantrunc.pcadim = 60;
matlabbatch{5}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{5}.spm.tools.beamforming.features.visualise = 1;
matlabbatch{6}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
matlabbatch{6}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
matlabbatch{7}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.woi = '<UNDEFINED>';
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.foi = '<UNDEFINED>';
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.contrast = [-1 1];
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.logpower = false;
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage';
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.scale = true;
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{7}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
matlabbatch{8}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{8}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
matlabbatch{8}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';
