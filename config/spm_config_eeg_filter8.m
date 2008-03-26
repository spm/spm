function S = spm_config_eeg_filter8
% configuration file for EEG Filtering
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_config_eeg_filter.m 1185 2008-03-04 16:31:21Z volkmar $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

typ = cfg_menu;
typ.tag = 'type';
typ.name = 'Filter type';
typ.labels = {'Butterworth'};
typ.values = {1};
typ.val = {1};
typ.help = {'Select the filter type.'};

PHz = cfg_entry;
PHz.tag = 'cutoff';
PHz.name = 'Cutoff';
PHz.strtype = 'r';
PHz.num = [1 1];
PHz.help = {'Enter the filter cutoff'};

flt = cfg_branch;
flt.tag = 'filter';
flt.name = 'Filter';
flt.val = {typ PHz};

S = cfg_exbranch;
S.tag = 'eeg_filter';
S.name = 'EEG Filter';
S.val = {D flt};
S.help = {'Low-pass filters EEG/MEG epoched data.'};
S.prog = @eeg_filter;
S.modality = {'EEG'};

function eeg_filter(job)
disp('hello')

% S.D = strvcat(S.D{:});
% spm_eeg_filter(S)
