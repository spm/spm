function convert = spm_cfg_eeg_convert
% Configuration file for M/EEG data conversion
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convert.m 6926 2016-11-09 22:13:19Z guillaume $

dataset = cfg_files;
dataset.tag = 'dataset';
dataset.name = 'File Name';
dataset.filter = 'any';
dataset.num = [1 1];
dataset.help = {'Select data set file.'};

timewin = cfg_entry;
timewin.tag = 'timewin';
timewin.name = 'Time window';
timewin.strtype = 'r';
timewin.num = [1 2];
timewin.help = {'start and end of epoch'};

readall = cfg_const;
readall.tag = 'readall';
readall.name = 'Read all';
readall.val  = {1};
readall.help = {''};

continuous = cfg_choice;
continuous.tag = 'continuous';
continuous.name = 'Continuous';
continuous.values = {timewin, readall};
continuous.val = {readall};
continuous.help = {''};

usetrials = cfg_const;
usetrials.tag = 'usetrials';
usetrials.name = 'Trials defined in data';
usetrials.val = {1};
usetrials.help = {''};

trlfile = cfg_files;
trlfile.tag = 'trlfile';
trlfile.name = 'Trial File';
trlfile.filter = 'mat';
trlfile.num = [1 1];
trlfile.help = {''};

conditionlabel = cfg_entry;
conditionlabel.tag = 'conditionlabel';
conditionlabel.name = 'Condition label';
conditionlabel.strtype = 's';
conditionlabel.help = {''};

eventtype = cfg_entry;
eventtype.tag = 'eventtype';
eventtype.name = 'Event type';
eventtype.strtype = 's';
eventtype.help = {''};

eventvalue = cfg_entry;
eventvalue.tag = 'eventvalue';
eventvalue.name = 'Event value';
eventvalue.strtype = 'r';
eventvalue.help = {''};

trialdef = cfg_branch;
trialdef.tag = 'trialdef';
trialdef.name = 'Trial';
trialdef.val = {conditionlabel, eventtype, eventvalue};
trialdef.help = {''};

define1 = cfg_repeat;
define1.tag = 'unused';
define1.name = 'Trial definitions';
define1.num = [1 Inf];
define1.values = {trialdef};
define1.help = {''};

define = cfg_branch;
define.tag = 'define';
define.name = 'Define trial';
define.val = {timewin define1};
define.help = {''};

epoched = cfg_choice;
epoched.tag = 'epoched';
epoched.name = 'Epoched';
epoched.values = {usetrials trlfile define};
epoched.help = {''};

mode = cfg_choice;
mode.tag = 'mode';
mode.name = 'Reading mode';
mode.values = {continuous, epoched};
mode.val = {continuous};
mode.help = {'Select whether you want to convert to continuous or epoched data.'};

outfile = cfg_entry;
outfile.tag = 'outfile';
outfile.name = 'Output filename';
outfile.strtype = 's';
outfile.num = [0 inf];
outfile.val = {''};
outfile.help = {'Choose filename. Leave empty to add ''spmeeg_'' to the input file name.'};

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
blocksize.help = {'size of blocks used internally to split large files default ~100Mb.'};

checkboundary = cfg_menu;
checkboundary.tag = 'checkboundary';
checkboundary.name = 'Check trial boundaries';
checkboundary.labels = {'Yes', 'No'};
checkboundary.val = {1};
checkboundary.values = {1,0};
checkboundary.help = {'Check if there are breaks in the file and do not read',...
    'across those breaks (recommended) or ignore them.'};

saveorigheader = cfg_menu;
saveorigheader.tag = 'saveorigheader';
saveorigheader.name = 'Save original header';
saveorigheader.labels = {'Yes', 'No'};
saveorigheader.val = {0};
saveorigheader.values = {1,0};
saveorigheader.help = {'Keep the original data header which might be useful for some later analyses.'};

inputformat = cfg_entry;
inputformat.tag = 'inputformat';
inputformat.name = 'Input data format';
inputformat.strtype = 's';
inputformat.val = {'autodetect'};
inputformat.num = [1 inf];
inputformat.help = {'Force the reader to assume a particular data format (usually not necessary).'};

convert = cfg_exbranch;
convert.tag = 'convert';
convert.name = 'Conversion';
convert.val = {dataset mode spm_cfg_eeg_channel_selector outfile...
    eventpadding blocksize checkboundary saveorigheader inputformat};
convert.help = {'Convert EEG/MEG data.'};
convert.prog = @eeg_convert;
convert.vout = @vout_eeg_convert;
convert.modality = {'EEG'};

function out = eeg_convert(job)

S = [];
S.dataset = job.dataset{1};
S.mode = job.mode;
S.channels = spm_cfg_eeg_channel_selector(job.channels);
if ~isempty(job.outfile)
    S.outfile = job.outfile;
end

S.eventpadding = job.eventpadding;
S.blocksize = job.blocksize;
S.checkboundary = job.checkboundary;
S.saveorigheader = job.saveorigheader;

if ~isequal(job.inputformat, 'autodetect')
    S.inputformat = job.inputformat;
end

S.mode = char((fieldnames(job.mode)));
switch  S.mode
    case 'continuous'
        if isfield(job.mode.continuous, 'timewin')
            S.timewin = job.mode.continuous.timewin;
        end
    case 'epoched'
        if isfield(job.mode.epoched, 'trlfile')
            S.trl = char(job.mode.epoched.trlfile);
        end
        
        if isfield(job.mode.epoched, 'define')
            S.trialdef =  job.mode.epoched.define.trialdef;
            S.timewin  =  job.mode.epoched.define.timewin;          
        end
end

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

