function S = spm_cfg_eeg_artefact
% configuration file for M/EEG artefact detection
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_artefact.m 3881 2010-05-07 21:02:57Z vladimir $

rev = '$Rev: 3881 $';

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
methods.val = {spm_cfg_eeg_channel_selector, fun};

methodsrep = cfg_repeat;
methodsrep.tag = 'methodsrep';
methodsrep.name = 'How to look for artefacts';
methodsrep.help = {'Choose channels and methods for artefact detection'};
methodsrep.values  = {methods};
methodsrep.num     = [1 Inf];

S = cfg_exbranch;
S.tag = 'artefact';
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
    S.methods(i).channels = spm_cfg_eeg_channel_selector(job.methods(i).channels);
    
    fun = fieldnames(job.methods(i).fun);
    fun = fun{1};
    
    S.methods(i).fun = fun;
    S.methods(i).settings = getfield(job.methods(i).fun, fun);
end
    
out.D = spm_eeg_artefact(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

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


