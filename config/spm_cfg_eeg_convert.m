function S = spm_cfg_eeg_filter
% configuration file for data conversion
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_convert.m 1562 2008-05-07 14:08:44Z stefan $

dataset = cfg_files;
dataset.tag = 'filename';
dataset.name = 'File Name';
dataset.filter = 'mat';
dataset.num = [1 1];
dataset.help = {'Select data set file.'};

continuous = cfg_menu;
continuous.tag = 'continuous';
continuous.name = 'Continuous data';
continuous.labels = {'Continuous', 'Epoched'};
continuous.values = {1, 0};
continuous.val = {1};
continuous.help = {'Select whether you want to convert to continuous or epoched data.'};

timewindow = cfg_entry;
timewindow.tag = 'timing';
timewindow.name = 'Timing';
timewindow.strtype = 'r';
timewindow.num = [1 2];
timewindow.help = {'start and end of epoch [s]'};

outfile = cfg_entry;
outfile.tag = 'Pout';
outfile.name = 'Output filename';
outfile.strtype = 's';
outfile.num = [1 inf];
outfile.help = {'Choose filename'};

% default for channels is for now 'all'

usetrials = cfg_menu;
usetrials.tag = 'usetrials';
usetrials.name = 'Use trials';
usetrials.labels = {'Use trials', 'Don''t use trials'};
usetrials.values = {1, 0};
usetrials.val = {1};
usetrials.help = {'1 - take the trials as defined in the data (default)',...
               '0 - use trial definition file even though the data is',...
                   'already epoched.'};
               
trlfile = cfg_entry;
trlfile.tag = 'trlfile';
trlfile.name = 'TRL filename';
trlfile.strtype = 's';
trlfile.num = [1 inf];
trlfile.help = {'Choose filename'};
     
datatype = cfg_menu;
datatype.tag = 'datatype';
datatype.name = 'Data type';
datatype.labels = {'float32','float64'};
datatype.val = {1};
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

conditionlabels = cfg_entry;
conditionlabels.tag = 'conditionlabels';
conditionlabels.name = 'Condition labels';
conditionlabels.strtype = 's';
conditionlabels.val = {'Undefined'};
conditionlabels.num = [1 inf];
conditionlabels.help = {'labels for the trials in the data Default - ''Undefined'''};
     
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
              
convert = cfg_branch;
convert.tag = 'convert';
convert.name = 'Convert';
convert.val = {dataset continuous timewindow outfile usetrials trlfile datatype eventpadding conditionlabels blocksize checkboundary};


S = cfg_exbranch;
S.tag = 'eeg_convert';
S.name = 'M/EEG Conversion';
S.val = {dataset convert};
S.help = {'Converts EEG/MEG data.'};
S.prog = @eeg_convert;
S.vout = @vout_eeg_convert;
S.modality = {'EEG'};

function out = eeg_convert(job)
% construct the S struct
S.D = job.D{1};
S.convert = job.convert;

out.D = spm_eeg_convert(S);

function dep = vout_eeg_convert(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Converted Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});
