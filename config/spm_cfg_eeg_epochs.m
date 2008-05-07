function S = spm_cfg_eeg_epochs
% configuration file for M/EEG epoching
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_epochs.m 1568 2008-05-07 18:23:23Z stefan $

rev = '$Rev: 1568 $';
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

gui = cfg_const;
gui.tag = 'gui';
gui.name = 'Use GUI';
gui.val{1} = double([]);

trl = cfg_entry;
trl.tag = 'trialdef';
trl.name = 'Trial matrix';
trl.strtype = 'r';
trl.num = [1 inf];
trl.help = {'Nx2 or Nx3 matrix (N - number of chantype) [start end offset]'};

cl = cfg_entry;
cl.tag = 'conditionlabels';
cl.name = 'Condition labels';
cl.strtype = 's';
cl.num = [inf inf];
cl.help = {'Enter condition labels'};

padding = cfg_entry;
padding.tag = 'padding';
padding.name = 'Padding';
padding.strtype = 'r';
padding.num = [1 1];
padding.help = {'Enter padding [s]'};

epochinfo = cfg_branch;
epochinfo.tag = 'epochinfo';
epochinfo.name = 'Epoch information';
epochinfo.val = {trl cl padding};

trialdef = cfg_entry;
trialdef.tag = 'trialdef';
trialdef.name = 'Trial definition';
trialdef.strtype = 'e';
trialdef.num = [inf inf];
trialdef.help = {'structure array for trial definition with fields',...
    'S.trialdef.conditionlabel - string label for the condition',...
    'S.trialdef.eventtype  - string',...
    'S.trialdef.eventvalue  - string, numeric or empty'};


trlchoice         = cfg_choice;
trlchoice.tag     = 'trialchoice';
trlchoice.name    = 'Choose a way how to define trials';
trlchoice.val = {gui};
trlchoice.help    = {'Choose one of the two options how to define trials'}';
trlchoice.values = {gui epochinfo trialdef};

events = cfg_branch;
events.tag = 'events';
events.name = 'Events';
events.val = {timing trlchoice};

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
S.pretrig = job.events.timing(1);
S.posttrig = job.events.timing(2);

if isfield(job.events.trlchoice, 'gui')

elseif isfield(job.events.trlchoice, 'trialdef')
    S.trialdef = 1;
    
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


