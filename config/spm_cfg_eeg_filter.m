function S = spm_cfg_eeg_filter
% configuration file for EEG Filtering
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_filter.m 4212 2011-02-23 17:50:55Z vladimir $

rev = '$Rev: 4212 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

typ = cfg_menu;
typ.tag = 'type';
typ.name = 'Filter type';
typ.labels = {'Butterworth', 'FIR'};
typ.values = {'butterworth', 'fir'};
typ.val = {'butterworth'};
typ.help = {'Select the filter type.'};

band = cfg_menu;
band.tag = 'band';
band.name = 'Filter band';
band.labels = {'Lowpass', 'Highpass', 'Bandpass', 'Stopband'};
band.values = {'low' 'high' 'bandpass' 'stop'};
band.val = {'low'};
band.help = {'Select the filter band.'};

PHz = cfg_entry;
PHz.tag = 'PHz';
PHz.name = 'Cutoff';
PHz.strtype = 'r';
PHz.num = [1 inf];
PHz.help = {'Enter the filter cutoff'};

dir = cfg_menu;
dir.tag = 'dir';
dir.name = 'Filter direction';
dir.labels = {'Zero phase', 'Forward', 'Backward'};
dir.values = {'twopass', 'onepass', 'onepass-reverse'};
dir.val = {'twopass'};
dir.help = {'Select the filter direction.'};

order = cfg_entry;
order.tag = 'order';
order.name = 'Filter order';
order.val = {5};
order.strtype = 'n';
order.num = [1 1];
order.help = {'Enter the filter order'};

flt = cfg_branch;
flt.tag = 'filter';
flt.name = 'Filter';
flt.val = {typ band PHz dir order};

S = cfg_exbranch;
S.tag = 'filter';
S.name = 'M/EEG Filter';
S.val = {D flt};
S.help = {'Low-pass filters EEG/MEG epoched data.'};
S.prog = @eeg_filter;
S.vout = @vout_eeg_filter;
S.modality = {'EEG'};

function out = eeg_filter(job)
% construct the S struct
S.D = job.D{1};
S.filter = job.filter;

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
