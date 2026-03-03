%-----------------------------------------------------------------------
% Job saved on 26-Aug-2025 13:43:11 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (00.00)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.convert.dataset = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
matlabbatch{1}.spm.meeg.convert.channels{1}.type = 'MEG';
matlabbatch{1}.spm.meeg.convert.channels{2}.type = 'EEG';
matlabbatch{1}.spm.meeg.convert.outfile = '';
matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Conversion: Converted Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.settype.channels{1}.chan = 'EEG061';
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.settype.channels{2}.chan = 'EEG062';
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.settype.newtype = 'EOG';
matlabbatch{2}.spm.meeg.preproc.prepare.task{2}.settype.channels{1}.chan = 'EEG063';
matlabbatch{2}.spm.meeg.preproc.prepare.task{2}.settype.newtype = 'ECG';
matlabbatch{2}.spm.meeg.preproc.prepare.task{3}.settype.channels{1}.chan = 'EEG064';
matlabbatch{2}.spm.meeg.preproc.prepare.task{3}.settype.newtype = 'Other';
matlabbatch{2}.spm.meeg.preproc.prepare.task{4}.setbadchan.channels{1}.chanfile = '<UNDEFINED>';
matlabbatch{2}.spm.meeg.preproc.prepare.task{4}.setbadchan.status = 1;
matlabbatch{3}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{3}.spm.meeg.preproc.downsample.fsample_new = 200;
matlabbatch{3}.spm.meeg.preproc.downsample.method = 'resample';
matlabbatch{3}.spm.meeg.preproc.downsample.prefix = 'd';
matlabbatch{4}.spm.meeg.preproc.filter.D(1) = cfg_dep('Downsampling: Downsampled Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{4}.spm.meeg.preproc.filter.type = 'butterworth';
matlabbatch{4}.spm.meeg.preproc.filter.band = 'high';
matlabbatch{4}.spm.meeg.preproc.filter.freq = 0.5;
matlabbatch{4}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{4}.spm.meeg.preproc.filter.order = 5;
matlabbatch{4}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{5}.spm.meeg.preproc.epoch.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{5}.spm.meeg.preproc.epoch.trialchoice.trlfile = '<UNDEFINED>';
matlabbatch{5}.spm.meeg.preproc.epoch.bc = 1;
matlabbatch{5}.spm.meeg.preproc.epoch.eventpadding = 0;
matlabbatch{5}.spm.meeg.preproc.epoch.prefix = 'e';
matlabbatch{6}.spm.meeg.preproc.bc.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{6}.spm.meeg.preproc.bc.timewin = [-100 0];
matlabbatch{6}.spm.meeg.preproc.bc.prefix = 'b';
matlabbatch{7}.spm.meeg.other.delete.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{8}.spm.meeg.other.delete.D(1) = cfg_dep('Downsampling: Downsampled Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{9}.spm.meeg.other.delete.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{10}.spm.meeg.other.delete.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
