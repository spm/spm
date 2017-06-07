function epoch = spm_cfg_eeg_epochs
% Configuration file for M/EEG epoching
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_epochs.m 7095 2017-06-07 10:15:59Z vladimir $


D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 1];
D.help   = {'Select the M/EEG mat file.'};

bc        = cfg_menu;
bc.tag    = 'bc';
bc.name   = 'Baseline correction';
bc.help   = {'Perform baseline correction when epoching, or not.'};
bc.labels = {'Yes','No'};
bc.values = {1 0};
bc.val    = {1};

% input via trl file
trlfile        = cfg_files;
trlfile.tag    = 'trlfile';
trlfile.name   = 'Trial definition file';
trlfile.filter = 'mat';
trlfile.num    = [1 1];
trlfile.help   = {'Select the trialfile mat file.'};

eventpadding         = cfg_entry;
eventpadding.tag     = 'eventpadding';
eventpadding.name    = 'Event padding';
eventpadding.strtype = 'r';
eventpadding.val     = {0};
eventpadding.num     = [1 1];
eventpadding.help    = {['In seconds: the additional time period around each trial',...
    'for which the events are saved with the trial (to let the ',...
    'user keep and use for analysis events which are outside ',...
    'trial borders). Default is 0 s.']};

% input via trialdef
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.help    = {'Start and end of epoch [ms].'};

conditionlabel         = cfg_entry;
conditionlabel.tag     = 'conditionlabel';
conditionlabel.name    = 'Condition label';
conditionlabel.strtype = 's';
conditionlabel.help    = {''};

eventtype         = cfg_entry;
eventtype.tag     = 'eventtype';
eventtype.name    = 'Event type';
eventtype.strtype = 's';
eventtype.help    = {''};

eventvalue         = cfg_entry;
eventvalue.tag     = 'eventvalue';
eventvalue.name    = 'Event value';
eventvalue.strtype = 'e';
eventvalue.help    = {''};

trlshift         = cfg_entry;
trlshift.tag     = 'trlshift';
trlshift.name    = 'Shift';
trlshift.strtype = 'r';
trlshift.num     = [1 1];
trlshift.val     = {0};
trlshift.help    = {'shift the triggers by a fixed amount (ms)',... 
                   '(e.g. projector delay).'};

trialdef      = cfg_branch;
trialdef.tag  = 'trialdef';
trialdef.name = 'Trial';
trialdef.val  = {conditionlabel eventtype eventvalue, trlshift};
trialdef.help = {''};

define1        = cfg_repeat;
define1.tag    = 'unused';
define1.name   = 'Trial definitions';
define1.values = {trialdef};
define1.help   = {''};

define      = cfg_branch;
define.tag  = 'define';
define.name = 'Define trial';
define.val  = {timewin define1};
define.help = {''};

% input via trialdef
trialength         = cfg_entry;
trialength.tag     = 'trialength';
trialength.name    = 'Trial length';
trialength.strtype = 'r';
trialength.num     = [1 1];
trialength.help    = {'Arbitary trial length [ms].'};

arbitrary      = cfg_branch;
arbitrary.tag  = 'arbitrary';
arbitrary.name = 'Arbitrary trials';
arbitrary.val  = {trialength, conditionlabel};
arbitrary.help = {'Epoch the data in arbitray fixed length trials (e.g. for spectral analysis).'};

trlchoice        = cfg_choice;
trlchoice.tag    = 'trialchoice';
trlchoice.name   = 'How to define trials';
trlchoice.help   = {'Choose one of the two options how to define trials.'}';
trlchoice.values = {trlfile define arbitrary};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the epoched dataset. Default prefix is ''e''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'e'};

epoch          = cfg_exbranch;
epoch.tag      = 'epoch';
epoch.name     = 'Epoching';
epoch.val      = {D, trlchoice, bc, eventpadding, prefix};
epoch.help     = {'Epoch continuous EEG/MEG data.'};
epoch.prog     = @eeg_epochs;
epoch.vout     = @vout_eeg_epochs;
epoch.modality = {'EEG'};


%==========================================================================
function out = eeg_epochs(job)
% construct the S struct
S.D = job.D{1};

if isfield(job.trialchoice, 'define')
    S.timewin  = job.trialchoice.define.timewin;
    S.trialdef = job.trialchoice.define.trialdef;
elseif isfield(job.trialchoice, 'arbitrary')
    S.trialength = job.trialchoice.arbitrary.trialength;
    S.conditionlabels = job.trialchoice.arbitrary.conditionlabel;
else    
    trlfile = load(char(job.trialchoice.trlfile));
    usetrl  = 0;
    
    % In the new code trl file contains both trl matrix and trial definition
    % struct. trl is usually only applicable to the file on which it was
    % defined and there it must have priority because it can be manually
    % adjusted in the GUI. Otherwise trialdef has priority if present.
    if isfield(trlfile, 'trl')
        if ~all(isfield(trlfile, {'trialdef', 'timewin'}))
             usetrl = 1;
        elseif isfield(trlfile, 'source')
                D = spm_eeg_load(S.D);
                if isequal(D.fname, trlfile.source)
                    usetrl = 1;
                end
        end
    elseif ~all(isfield(trlfile, {'trialdef', 'timewin'}))
        error('The trial definition file could not be inetrpreted');
    end
    
    if usetrl
        S.trl = trlfile.trl;
        if isfield(trlfile, 'conditionlabels')
            S.conditionlabels = trlfile.conditionlabels;
        end
    else
        S.trialdef = trlfile.trialdef;
        S.timewin = trlfile.timewin;
    end
end

S.bc = job.bc;
S.prefix = job.prefix;
S.eventpadding = job.eventpadding;

out.D      = spm_eeg_epochs(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};


%==========================================================================
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
