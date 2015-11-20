function downsample = spm_cfg_eeg_downsample
% Configuration file for M/EEG downsampling
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_downsample.m 6602 2015-11-20 19:04:49Z vladimir $


D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 1];
D.help   = {'Select the M/EEG mat file.'};

fsample_new         = cfg_entry;
fsample_new.tag     = 'fsample_new';
fsample_new.name    = 'New sampling rate';
fsample_new.strtype = 'r';
fsample_new.num     = [1 1];
fsample_new.help    = {'Input the new sampling rate [Hz].'};

method = cfg_menu;
method.tag  = 'method';
method.name = 'Resampling method';
method.labels = {'resample', 'decimate', 'downsample', 'fft'};
method.values = {'resample', 'decimate', 'downsample', 'fft'};
method.val = {'resample'};
method.help = {'Select the downsampling method.'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the filtered dataset. Default prefix is ''d''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'d'};

downsample      = cfg_exbranch;
downsample.tag  = 'downsample';
downsample.name = 'Downsampling';
downsample.val  = {D fsample_new, method, prefix};
downsample.help = {'Downsample EEG/MEG data.'};
downsample.prog = @eeg_downsample;
downsample.vout = @vout_eeg_downsample;
downsample.modality = {'EEG'};


%==========================================================================
function out = eeg_downsample(job)
S   = job;
S.D = S.D{1};

out.D = spm_eeg_downsample(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};


%==========================================================================
function dep = vout_eeg_downsample(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Downsampled data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Downsampled Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
