function S = spm_cfg_eeg_epochs
% configuration file for M/EEG epoching
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_epochs.m 3881 2010-05-07 21:02:57Z vladimir $

rev = '$Rev: 3881 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

% input via trl file
trlfile = cfg_files;
trlfile.tag = 'trlfile';
trlfile.name = 'File Name';
trlfile.filter = 'mat';
trlfile.num = [1 1];
trlfile.help = {'Select the trialfile mat file.'};

padding = cfg_entry;
padding.tag = 'padding';
padding.name = 'Padding';
padding.strtype = 'r';
padding.num = [1 1];
padding.help = {'Enter padding [s]: the additional time period around each trial',...
    'for which the events are saved with the trial (to let the',...
    'user keep and use for analysis events which are outside'};

epochinfo = cfg_branch;
epochinfo.tag = 'epochinfo';
epochinfo.name = 'Epoch information';
epochinfo.val = {trlfile padding};

% input via trialdef
timewindow = cfg_entry;
timewindow.tag = 'timewindow';
timewindow.name = 'Timing';
timewindow.strtype = 'r';
timewindow.num = [1 2];
timewindow.help = {'start and end of epoch [ms]'};

conditionlabel = cfg_entry;
conditionlabel.tag = 'conditionlabel';
conditionlabel.name = 'Condition label';
conditionlabel.strtype = 's';

eventtype = cfg_entry;
eventtype.tag = 'eventtype';
eventtype.name = 'Event type';
eventtype.strtype = 's';

eventvalue = cfg_entry;
eventvalue.tag = 'eventvalue';
eventvalue.name = 'Event value';
eventvalue.strtype = 'e';

trialdef = cfg_branch;
trialdef.tag = 'trialdef';
trialdef.name = 'Trial';
trialdef.val = {conditionlabel, eventtype, eventvalue};

define1 = cfg_repeat;
define1.tag = 'unused';
define1.name = 'Trial definitions';
define1.values = {trialdef};

define = cfg_branch;
define.tag = 'define';
define.name = 'Define trial';
define.val = {timewindow define1};

trlchoice         = cfg_choice;
trlchoice.tag     = 'trialchoice';
trlchoice.name    = 'Choose a way how to define trials';
trlchoice.help    = {'Choose one of the two options how to define trials'}';
trlchoice.values = {epochinfo define};


S = cfg_exbranch;
S.tag = 'epoch';
S.name = 'M/EEG Epoching';
S.val = {D trlchoice};
S.help = {'Epoch continuous EEG/MEG data.'};
S.prog = @eeg_epochs;
S.vout = @vout_eeg_epochs;
S.modality = {'EEG'};


function out = eeg_epochs(job)
% construct the S struct
S.D = job.D{1};

if isfield(job.trialchoice, 'define')
    S.pretrig = job.trialchoice.define.timewindow(1);
    S.posttrig = job.trialchoice.define.timewindow(2);
    S.trialdef = job.trialchoice.define.trialdef;
else
    S.epochinfo = job.trialchoice.epochinfo;
    S.epochinfo.trlfile = S.epochinfo.trlfile{1};
end

% set review and save options both to 0 to not pop up something
S.reviewtrials = 0;
S.save = 0;

out.D = spm_eeg_epochs(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_epochs(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Epoched Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Epoched Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

