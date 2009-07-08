function S = spm_cfg_eeg_artefact
% configuration file for M/EEG artefact detection
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_artefact.m 3258 2009-07-08 17:46:54Z vladimir $

rev = '$Rev: 3258 $';

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

badchanthresh = cfg_entry;
badchanthresh.tag = 'badchanthresh';
badchanthresh.name = 'Bad channel threshold';
badchanthresh.strtype = 'r';
badchanthresh.num = [1 1];
badchanthresh.val = {0.2};
badchanthresh.help = {'Fraction of trials with artefacts ', ...
    'above which an M/EEG channel is declared as bad.'};

chanall = cfg_const;
chanall.tag = 'type';
chanall.name = 'All';
chanall.val = {'all'};

chanmeg = cfg_const;
chanmeg.tag = 'type';
chanmeg.name = 'MEG';
chanmeg.val = {'MEG'};

chanmegplanar = cfg_const;
chanmegplanar.tag = 'type';
chanmegplanar.name = 'MEGPLANAR';
chanmegplanar.val = {'MEGPLANAR'};

chaneeg = cfg_const;
chaneeg.tag = 'type';
chaneeg.name = 'EEG';
chaneeg.val = {'EEG'};

chaneog = cfg_const;
chaneog.tag = 'type';
chaneog.name = 'EOG';
chaneog.val = {'EOG'};

chanecg = cfg_const;
chanecg.tag = 'type';
chanecg.name = 'ECG';
chanecg.val = {'ECG'};

chanemg = cfg_const;
chanemg.tag = 'type';
chanemg.name = 'EMG';
chanemg.val = {'EMG'};

chanlfp = cfg_const;
chanlfp.tag = 'type';
chanlfp.name = 'LFP';
chanlfp.val = {'LFP'};

chanfile = cfg_files;
chanfile.tag = 'file';
chanfile.name = 'Channel file';
chanfile.filter = 'mat';
chanfile.num = [1 1];

channels = cfg_choice;
channels.tag = 'channels';
channels.name = 'Channel selection';
channels.values = {chanall, chanmeg, chanmegplanar, chaneeg, chaneog, chanecg, chanemg, chanlfp, chanfile};
channels.val = {chanall};

artefact_funs = dir(fullfile(spm('dir'), 'spm_eeg_artefact_*.m'));
artefact_funs = {artefact_funs(:).name};

fun      = cfg_choice;
fun.tag  = 'fun';
fun.name = 'Detection algorithm';
for i = 1:numel(artefact_funs)
    fun.values{i} = feval(spm_str_manip(artefact_funs{i}, 'r'));
end

methods = cfg_branch;
methods.tag = 'methods';
methods.name = 'Method';
methods.val = {channels, fun};

methodsrep = cfg_repeat;
methodsrep.tag = 'methodsrep';
methodsrep.name = 'How to look for artefacts';
methodsrep.help = {'Choose channels and methods for artefact detection'};
methodsrep.values  = {methods};
methodsrep.num     = [1 Inf];


S = cfg_exbranch;
S.tag = 'eeg_artefact';
S.name = 'M/EEG Artefact detection';
S.val = {D, badchanthresh, methodsrep};
S.help = {'Detect artefacts in epoched M/EEG data.'};
S.prog = @eeg_artefact;
S.vout = @vout_eeg_artefact;
S.modality = {'EEG'};


function out = eeg_artefact(job)
% construct the S struct
S.D = job.D{1};
S.badchanthresh = job.badchanthresh;

for i = 1:numel(job.methods)
    if isfield(job.methods(i).channels, 'type')
        S.methods(i).channels = job.methods(i).channels.type;
    else
        S.methods(i).channels = getfield(load(job.methods(i).channels.file{1}), 'label');
    end
    
    fun = fieldnames(job.methods(i).fun);
    fun = fun{1};
    
    S.methods(i).fun = fun;
    S.methods(i).settings = getfield(job.methods(i).fun, fun);
end
    
out.D = spm_eeg_artefact(S);
out.Dfname = {out.D.fname};

function dep = vout_eeg_artefact(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Artefact detection';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Artefact-detected Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


