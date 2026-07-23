% DAiSS HEP example pipeline: LCMV source waveform reconstruction for the 
% heartbeat-evoked potential (HEP) in EEG data.
%
% Extracts virtual-electrode time courses from two VOIs using an LCMV
% beamformer, writes them back into an SPM M/EEG dataset, then averages
% over trials. 
%
% Example VOIs in step 3 (to be modified for your own VOIs):
%   Insula     [36  0  0]
%   Right ACC  [ 6 40 10]
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

% Step 3: Define sources (VOIs)
matlabbatch{3}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{3}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{1}.voidef.label = 'insula';
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{1}.voidef.pos = [36 0 0];
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{1}.voidef.ori = [0 0 0];
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{2}.voidef.label = 'rightACC';
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{2}.voidef.pos = [6 40 10];
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.vois{2}.voidef.ori = [0 0 0];
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.radius = 20;
matlabbatch{3}.spm.tools.beamforming.sources.plugin.voi.resolution = 5;
matlabbatch{3}.spm.tools.beamforming.sources.normalise_lf = false;
matlabbatch{3}.spm.tools.beamforming.sources.visualise = 0;

% Step 4: Covariance 
matlabbatch{4}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{4}.spm.tools.beamforming.features.woi = [-Inf Inf];
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

% Step 6: Output montage
matlabbatch{6}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.output.plugin.montage.method = 'max';
matlabbatch{6}.spm.tools.beamforming.output.plugin.montage.vois = cell(1, 0);

% Step 7: Write source waveforms to SPM M/EEG file
matlabbatch{7}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','BF'));
matlabbatch{7}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'onlinecopy';
matlabbatch{7}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'EEG';
matlabbatch{7}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
matlabbatch{7}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';

% Step 8: Average
matlabbatch{8}.spm.meeg.averaging.average.D(1) = cfg_dep('Write: Output files', ...
    substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{8}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{8}.spm.meeg.averaging.average.plv = false;
matlabbatch{8}.spm.meeg.averaging.average.prefix = 'm';