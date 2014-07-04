%-----------------------------------------------------------------------
% Job saved on 04-Jul-2014 13:55:33 by cfg_util (rev $Rev: 6090 $)
% spm SPM - SPM12b (beta)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.convert.dataset = {'subject1.bdf'};
matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
matlabbatch{1}.spm.meeg.convert.channels{1}.type = 'EEG';
matlabbatch{1}.spm.meeg.convert.outfile = '';
matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Conversion: Converted Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.avref.fname = 'avref_montage.mat';
matlabbatch{3}.spm.meeg.preproc.montage.D(1) = cfg_dep('Conversion: Converted Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.montspec.montage.montagefile(1) = cfg_dep('Prepare: Average reference montage', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avrefname'));
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.montspec.montage.keepothers = 0;
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.blocksize = 655360;
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.prefix = 'M';
matlabbatch{4}.spm.meeg.preproc.filter.D(1) = cfg_dep('Montage: Montaged Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{4}.spm.meeg.preproc.filter.type = 'butterworth';
matlabbatch{4}.spm.meeg.preproc.filter.band = 'high';
matlabbatch{4}.spm.meeg.preproc.filter.freq = 0.1;
matlabbatch{4}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{4}.spm.meeg.preproc.filter.order = 5;
matlabbatch{4}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{5}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{5}.spm.meeg.preproc.downsample.fsample_new = 200;
matlabbatch{5}.spm.meeg.preproc.downsample.prefix = 'd';
matlabbatch{6}.spm.meeg.preproc.filter.D(1) = cfg_dep('Downsampling: Downsampled Datafile', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{6}.spm.meeg.preproc.filter.type = 'butterworth';
matlabbatch{6}.spm.meeg.preproc.filter.band = 'low';
matlabbatch{6}.spm.meeg.preproc.filter.freq = 30;
matlabbatch{6}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{6}.spm.meeg.preproc.filter.order = 5;
matlabbatch{6}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{7}.spm.meeg.preproc.epoch.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-100 400];
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'std';
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 1;
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'odd';
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'STATUS';
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 3;
matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
matlabbatch{7}.spm.meeg.preproc.epoch.bc = 1;
matlabbatch{7}.spm.meeg.preproc.epoch.eventpadding = 0;
matlabbatch{7}.spm.meeg.preproc.epoch.prefix = 'e';
matlabbatch{8}.spm.meeg.preproc.artefact.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{8}.spm.meeg.preproc.artefact.mode = 'reject';
matlabbatch{8}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
matlabbatch{8}.spm.meeg.preproc.artefact.append = true;
matlabbatch{8}.spm.meeg.preproc.artefact.methods.channels{1}.type = 'EEG';
matlabbatch{8}.spm.meeg.preproc.artefact.methods.fun.threshchan.threshold = 80;
matlabbatch{8}.spm.meeg.preproc.artefact.methods.fun.threshchan.excwin = 1000;
matlabbatch{8}.spm.meeg.preproc.artefact.prefix = 'a';
matlabbatch{9}.spm.meeg.averaging.average.D(1) = cfg_dep('Artefact detection: Artefact-detected Datafile', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{9}.spm.meeg.averaging.average.userobust.robust.ks = 3;
matlabbatch{9}.spm.meeg.averaging.average.userobust.robust.bycondition = false;
matlabbatch{9}.spm.meeg.averaging.average.userobust.robust.savew = false;
matlabbatch{9}.spm.meeg.averaging.average.userobust.robust.removebad = true;
matlabbatch{9}.spm.meeg.averaging.average.plv = false;
matlabbatch{9}.spm.meeg.averaging.average.prefix = 'm';

spm_jobman('run', matlabbatch);
