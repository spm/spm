%-----------------------------------------------------------------------
% Job saved on 28-Jan-2026 16:18:07 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (00.00)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.preproc.crop.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.preproc.crop.timewin = [-100 800];
matlabbatch{1}.spm.meeg.preproc.crop.freqwin = [-Inf Inf];
matlabbatch{1}.spm.meeg.preproc.crop.channels{1}.all = 'all';
matlabbatch{1}.spm.meeg.preproc.crop.prefix = 'p';
matlabbatch{2}.spm.meeg.preproc.artefact.D(1) = cfg_dep('Crop: Cropped M/EEG datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.artefact.mode = 'reject';
matlabbatch{2}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
matlabbatch{2}.spm.meeg.preproc.artefact.append = true;
matlabbatch{2}.spm.meeg.preproc.artefact.methods.channels{1}.type = 'EOG';
matlabbatch{2}.spm.meeg.preproc.artefact.methods.fun.threshchan.threshold = 200;
matlabbatch{2}.spm.meeg.preproc.artefact.methods.fun.threshchan.excwin = 1000;
matlabbatch{2}.spm.meeg.preproc.artefact.prefix = 'a';
matlabbatch{3}.spm.meeg.averaging.average.D(1) = cfg_dep('Artefact detection: Artefact-detected Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{3}.spm.meeg.averaging.average.userobust.standard = false;
matlabbatch{3}.spm.meeg.averaging.average.plv = false;
matlabbatch{3}.spm.meeg.averaging.average.prefix = 'm';
matlabbatch{4}.spm.meeg.averaging.contrast.D(1) = cfg_dep('Averaging: Averaged Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(1).c = [0.5 0.5 -1];
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(1).label = 'Faces_vs_Scrambled';
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(2).c = [1 -1 0];
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(2).label = 'Famous_vs_Unfamiliar';
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(3).c = [1 0 0];
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(3).label = 'Famous';
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(4).c = [0 1 0];
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(4).label = 'Unfamiliar';
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(5).c = [0 0 1];
matlabbatch{4}.spm.meeg.averaging.contrast.contrast(5).label = 'Scrambled';
matlabbatch{4}.spm.meeg.averaging.contrast.weighted = 0;
matlabbatch{4}.spm.meeg.averaging.contrast.prefix = 'w';
matlabbatch{5}.spm.meeg.other.delete.D(1) = cfg_dep('Averaging: Averaged Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{6}.spm.meeg.images.convert2images.D(1) = cfg_dep('Contrast over epochs: Contrast Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{6}.spm.meeg.images.convert2images.mode = 'scalp x time';
matlabbatch{6}.spm.meeg.images.convert2images.conditions = {};
matlabbatch{6}.spm.meeg.images.convert2images.channels{1}.type = 'EEG';
matlabbatch{6}.spm.meeg.images.convert2images.optimise = 1;
matlabbatch{6}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
matlabbatch{6}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
matlabbatch{6}.spm.meeg.images.convert2images.prefix = 'eeg_evoked_';
matlabbatch{7}.spm.meeg.images.convert2images.D(1) = cfg_dep('Contrast over epochs: Contrast Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{7}.spm.meeg.images.convert2images.mode = 'scalp x time';
matlabbatch{7}.spm.meeg.images.convert2images.conditions = {};
matlabbatch{7}.spm.meeg.images.convert2images.channels{1}.type = 'MEG';
matlabbatch{7}.spm.meeg.images.convert2images.optimise = 1;
matlabbatch{7}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
matlabbatch{7}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
matlabbatch{7}.spm.meeg.images.convert2images.prefix = 'meg_evoked_';
