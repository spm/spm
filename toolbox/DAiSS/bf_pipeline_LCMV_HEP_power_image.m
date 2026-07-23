% DAiSS example pipeline: LCMV beamforming source power-contrast imaging 
% for the heartbeat-evoked potential (HEP) in EEG data.
%
% Computes source power in a baseline and an activation time window, 
% writes each as a NIfTI, and forms the activation-minus-baseline 
% difference image.
%
% Example time windows (to be modified for your own windows of interest):
%   Covariance : [-100 400] ms (this should include the QRS window 
% and HEP window)
%   Baseline   : [-100 100] ms (this should be the QRS window)
%   Activation : [ 200 400] ms (this shoudl be the HEP window)
%
% Note: this code was validated on simulated EEG data with a known 
% HEP window. In empricial recordings, source-level SNR is typically too 
% low to yield meaningful results.
%
% To use: select this pipeline from the DAiSS menu (spm_beamforming),
% then set the input dataset (module 1, "M/EEG dataset") and the output
% directory (module 2, "Directory") in the Batch Editor before running.

% To find out more about this pipeline and the validation results, 
% please refer to the following preprint:
% Virjee, R.-I., Kandasamy, R., Garfinkel, S.N., Yogarajah, M., Litvak, V. 
% and Carmichael, D.W. (2026). 
% Can the heartbeat-evoked potential (HEP) be separated from Cardiac 
% Artefact (CA) using beamforming? doi:10.64898/2026.07.11.737958.

% Vladimir Litvak and Rania-Iman Virjee
% Copyright (C) 2026 Department of Imaging Neuroscience


% Step 1: Head model
matlabbatch{1}.spm.meeg.source.headmodel.D = {''};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregdefault = 1;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

% Step 2: Prepare data
matlabbatch{2}.spm.tools.beamforming.data.dir = '<UNDEFINED>';
matlabbatch{2}.spm.tools.beamforming.data.D(1) = cfg_dep('M/EEG head model specification: M/EEG dataset(s) with a forward model', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','D'));
matlabbatch{2}.spm.tools.beamforming.data.val = 1;
matlabbatch{2}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{2}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{2}.spm.tools.beamforming.data.overwrite = 0;

% Step 3: Define sources
matlabbatch{3}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{3}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.grid.resolution = 5;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';
matlabbatch{3}.spm.tools.beamforming.sources.plugin.grid.constrain = 'iskull';
matlabbatch{3}.spm.tools.beamforming.sources.normalise_lf = false;
matlabbatch{3}.spm.tools.beamforming.sources.visualise = 0;

% Step 4: Covariance
matlabbatch{4}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{4}.spm.tools.beamforming.features.woi = [-100 400];
matlabbatch{4}.spm.tools.beamforming.features.modality = {'EEG'};
matlabbatch{4}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{4}.spm.tools.beamforming.features.cross_terms = 'megeeg';
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.foi = [0 30];
matlabbatch{4}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
matlabbatch{4}.spm.tools.beamforming.features.regularisation.mantrunc.pcadim = 100;
matlabbatch{4}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{4}.spm.tools.beamforming.features.visualise = 0;

% Step 5: Inverse 
matlabbatch{5}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
matlabbatch{5}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;

% Step 6: Baseline-window power image
matlabbatch{6}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.woi = [-100 100];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.foi = [0 30];
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.erp = true;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.logpower = false;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.scale = true;
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{6}.spm.tools.beamforming.output.plugin.image_power.modality = 'EEG';

% Step 7: Write baseline image
matlabbatch{7}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{7}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
matlabbatch{7}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';

% Step 8: Move baseline
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Write: Output files', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = {pwd};
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = '.*';
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'baseline';
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = true;

% Step 9: Activation-window power image
matlabbatch{9}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.woi = [200 400];
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.foi = [0 30];
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.erp = true;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.logpower = false;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage';
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.scale = true;
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{9}.spm.tools.beamforming.output.plugin.image_power.modality = 'EEG';

% Step 10: Write activation
matlabbatch{10}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', ...
    substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{10}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
matlabbatch{10}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';

% Step 11: Move activation
matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Write: Output files', ...
    substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = {pwd};
matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = '.*';
matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'activation';
matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = true;

% Step 12: Calculate difference
matlabbatch{12}.spm.util.imcalc.input(1) = cfg_dep('Move/Delete Files: Moved/Copied Files', ...
    substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{12}.spm.util.imcalc.input(2) = cfg_dep('Move/Delete Files: Moved/Copied Files', ...
    substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{12}.spm.util.imcalc.output = 'output';
matlabbatch{12}.spm.util.imcalc.outdir = {''};
matlabbatch{12}.spm.util.imcalc.expression = 'i2-i1';
matlabbatch{12}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{12}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{12}.spm.util.imcalc.options.mask = 0;
matlabbatch{12}.spm.util.imcalc.options.interp = 1;
matlabbatch{12}.spm.util.imcalc.options.dtype = 16;

