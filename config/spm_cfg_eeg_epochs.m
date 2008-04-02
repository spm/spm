function S = spm_cfg_eeg_epochs
% configuration file for M/EEG epoching
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_epochs.m 1295 2008-04-02 14:31:24Z volkmar $

rev = '$Rev: 1295 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

timing = cfg_entry;
timing.tag = 'timing';
timing.name = 'Timing';
timing.strtype = 'r';
timing.num = [1 2];
timing.help = {'start and end of epoch [ms]'};

types = cfg_entry;
types.tag = 'types';
types.name = 'Trial types';
types.strtype = 'n';
types.num = [1 inf];
types.help = {'Select the indices of trial types to be epoched.'};

Ec = cfg_entry;
Ec.tag = 'newcodes';
Ec.name = 'New event codes';
Ec.strtype = 'r';
Ec.num = [0 inf];
Ec.help = {'New event codes'};

nothing = cfg_const;
nothing.tag = 'nothing';
nothing.name = 'No new event codes';
nothing.val{1} = double([]);

Inewlist         = cfg_choice;
Inewlist.tag     = 'Inewlist';
Inewlist.name    = 'Use new list of events';
Inewlist.val = {nothing};
Inewlist.help    = {'Choose this option if you want to assign new'
                     'event codes.'}';
Inewlist.values = {Ec nothing};

newlabels = cfg_entry;
newlabels.tag = 'newlabels';
newlabels.name = 'New labels';
newlabels.strtype = 'e';
newlabels.num = [inf inf];
newlabels.help = {
    'Add one label (a string) for each epoched'
    'trial type'};

events = cfg_branch;
events.tag = 'events';
events.name = 'Events';
events.val = {timing types Inewlist newlabels};

S = cfg_exbranch;
S.tag = 'eeg_epochs';
S.name = 'M/EEG epoching';
S.val = {D events};
S.help = {'Epoch continuous EEG/MEG data.'};
S.prog = @eeg_epochs;
S.vout = @vout_eeg_epochs;
S.modality = {'EEG'};


function out = eeg_epochs(job)
% construct the S struct
S.D = job.D{1};
S.events = job.events;
S.events.start = S.events.timing(1);
S.events.stop = S.events.timing(2);
if isfield(S.events.Inewlist, 'nothing')
    S.events.Inewlist = 0;
else
    S.events.Inewlist = 1;
    S.events.Ec = job.events.Inewlist.newcodes;
end

out.D = spm_eeg_epochs(S);

function dep = vout_eeg_epochs(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Epoched Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


