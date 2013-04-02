function filter = spm_cfg_eeg_filter
% configuration file for EEG Filtering
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_filter.m 5377 2013-04-02 17:07:57Z vladimir $

rev = '$Rev: 5377 $';

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

type = cfg_menu;
type.tag = 'type';
type.name = 'Type';
type.labels = {'Butterworth', 'FIR'};
type.values = {'butterworth', 'fir'};
type.val = {'butterworth'};
type.help = {'Select the filter typee.'};

band = cfg_menu;
band.tag = 'band';
band.name = 'Band';
band.labels = {'Lowpass', 'Highpass', 'Bandpass', 'Stopband'};
band.values = {'low' 'high' 'bandpass' 'stop'};
band.val = {'low'};
band.help = {'Select the filter band.'};

freq = cfg_entry;
freq.tag = 'freq';
freq.name = 'Cutoff(s)';
freq.strtype = 'r';
freq.num = [1 inf];
freq.help = {'Enter the filter cutoff'};

dir = cfg_menu;
dir.tag = 'dir';
dir.name = 'Direction';
dir.labels = {'Zero phase', 'Forward', 'Backward'};
dir.values = {'twopass', 'onepass', 'onepass-reverse'};
dir.val = {'twopass'};
dir.help = {'Select the filter direction.'};

order = cfg_entry;
order.tag = 'order';
order.name = 'Order';
order.val = {5};
order.strtype = 'n';
order.num = [1 1];
order.help = {'Enter the filter order'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the filtered dataset. Default prefix is ''f''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'f'};

filter = cfg_exbranch;
filter.tag = 'filter';
filter.name = 'Filter';
filter.val = {D type band freq dir order, prefix};
filter.help = {'Filters EEG/MEG data.'};
filter.prog = @eeg_filter;
filter.vout = @vout_eeg_filter;
filter.modality = {'EEG'};

function out = eeg_filter(job)
% construct the S struct
S = job;
S.D = S.D{1};

out.D = spm_eeg_filter(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_filter(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Filtered Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Filtered Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
