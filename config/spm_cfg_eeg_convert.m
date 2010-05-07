function S = spm_cfg_eeg_convert
% configuration file for data conversion
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convert.m 3881 2010-05-07 21:02:57Z vladimir $

dataset = cfg_files;
dataset.tag = 'dataset';
dataset.name = 'File Name';
dataset.filter = 'any';
dataset.num = [1 1];
dataset.help = {'Select data set file.'};

timewindow = cfg_entry;
timewindow.tag = 'timing';
timewindow.name = 'Timing';
timewindow.strtype = 'r';
timewindow.num = [1 2];
timewindow.help = {'start and end of epoch [ms]'};

readall = cfg_const;
readall.tag = 'readall';
readall.name = 'Read all';
readall.val  = {1};

read = cfg_choice;
read.tag = 'read';
read.name = 'Continuous';
read.values = {timewindow, readall};
read.val = {readall};

usetrials = cfg_const;
usetrials.tag = 'usetrials';
usetrials.name = 'Trials defined in data';
usetrials.val = {1};

trlfile = cfg_files;
trlfile.tag = 'trlfile';
trlfile.name = 'Trial File';
trlfile.filter = 'mat';
trlfile.num = [1 1];

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

trials = cfg_choice;
trials.tag = 'trials';
trials.name = 'Epoched';
trials.values = {usetrials trlfile define};

continuous = cfg_choice;
continuous.tag = 'continuous';
continuous.name = 'Reading mode';
continuous.values = {read trials};
continuous.val = {read};
continuous.help = {'Select whether you want to convert to continuous or epoched data.'};

chanall = cfg_const;
chanall.tag = 'chanall';
chanall.name = 'All';
chanall.val = {1};

chanmeg = cfg_const;
chanmeg.tag = 'chanmeg';
chanmeg.name = 'MEG';
chanmeg.val = {1};

chaneeg = cfg_const;
chaneeg.tag = 'chaneeg';
chaneeg.name = 'EEG';
chaneeg.val = {1};

chanfile = cfg_files;
chanfile.tag = 'chanfile';
chanfile.name = 'Channel file';
chanfile.filter = 'mat';
chanfile.num = [1 1];

channels = cfg_choice;
channels.tag = 'channels';
channels.name = 'Channel selection';
channels.values = {chanall,chanmeg,chaneeg,chanfile};
channels.val = {chanall};

outfile = cfg_entry;
outfile.tag = 'outfile';
outfile.name = 'Output filename';
outfile.strtype = 's';
outfile.num = [0 inf];
outfile.val = {''};
outfile.help = {'Choose filename. Leave empty to add ''spm8_'' to the input file name.'};

datatype = cfg_menu;
datatype.tag = 'datatype';
datatype.name = 'Data type';
datatype.labels = {'float32-le','float64-le'};
datatype.val    = {'float32-le'};
datatype.values = {'float32-le','float64-le'};
datatype.help = {'Determine data type to save data in.'};

eventpadding = cfg_entry;
eventpadding.tag = 'eventpadding';
eventpadding.name = 'Event padding';
eventpadding.strtype = 'r';
eventpadding.val = {0};
eventpadding.num = [1 1];
eventpadding.help = {'in sec - the additional time period around each trial',...
    'for which the events are saved with the trial (to let the',...
    'user keep and use for analysis events which are outside',...
    'trial borders). Default - 0'};

blocksize = cfg_entry;
blocksize.tag = 'blocksize';
blocksize.name = 'Block size';
blocksize.strtype = 'r';
blocksize.val = {3276800};
blocksize.num = [1 1];
blocksize.help = {'size of blocks used internally to split large files default ~100Mb'};

checkboundary = cfg_menu;
checkboundary.tag = 'checkboundary';
checkboundary.name = 'Check boundary';
checkboundary.labels = {'Check boundaries', 'Don''t check boundaries'};
checkboundary.val = {1};
checkboundary.values = {1,0};
checkboundary.help = {'1 - check if there are breaks in the file and do not read',...
    'across those breaks (default).',...
    '0 - ignore breaks (not recommended)'};

inputformat = cfg_entry;
inputformat.tag = 'inputformat';
inputformat.name = 'Input data format';
inputformat.strtype = 's';
inputformat.val = {'autodetect'};
inputformat.num = [1 inf];
inputformat.help = {'Force the reader to assume a particular data format (usually not necessary)'};

S = cfg_exbranch;
S.tag = 'convert';
S.name = 'M/EEG Conversion';
S.val = {dataset continuous channels outfile datatype eventpadding blocksize checkboundary inputformat};
S.help = {'Converts EEG/MEG data.'};
S.prog = @eeg_convert;
S.vout = @vout_eeg_convert;
S.modality = {'EEG'};

function out = eeg_convert(job)
% construct the S struct
S.dataset = job.dataset{1};
S.continuous = job.continuous;
S.channels = job.channels;
if ~isempty(job.outfile)
    S.outfile = job.outfile;
end
S.datatype = job.datatype;
S.eventpadding = job.eventpadding;
S.blocksize = job.blocksize;
S.checkboundary = job.checkboundary;

if ~isequal(job.inputformat, 'autodetect')
    S.inputformat = job.inputformat;
end

if isfield(S.continuous, 'read')
    S.continuous = 1;
    if isfield(job.continuous.read, 'timing')
        S.timewindow = job.continuous.read.timing;
    end
else
    if isfield(S.continuous.trials, 'usetrials')
        S.usetrials = S.continuous.trials.usetrials;
    end
    if isfield(S.continuous.trials, 'trlfile')
        S.trlfile = char(S.continuous.trials.trlfile);
        S.usetrials = 0;
    end
    
    if isfield(S.continuous.trials, 'define')
        S.trialdef = S.continuous.trials.define.trialdef;
        S.pretrig = S.continuous.trials.define.timing(1);
        S.posttrig = S.continuous.trials.define.timing(2);
        S.reviewtrials = 0;
        S.save = 0;
        S.usetrials = 0;
        [S.trl, S.conditionlabel] = spm_eeg_definetrial(S);
    end
    S.continuous = 0;
end

if isfield(S.channels, 'chanmeg')
    S.channels = 'MEG';
elseif isfield(S.channels, 'chaneeg')
    S.channels = 'EEG';
elseif isfield(S.channels, 'chanall')
    S.channels = 'all';
elseif isfield(S.channels, 'chanfile')
    S.chanfile = S.channels.chanfile{1};
    S.channels = 'file';
end
S.save = 0;
S.review = 0;

S = spm_eeg_channelselection(S);
S = rmfield(S, 'save');
S = rmfield(S, 'review');

out.D = spm_eeg_convert(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_convert(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Converted M/EEG Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Converted Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

