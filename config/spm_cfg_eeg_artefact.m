function artefact = spm_cfg_eeg_artefact
% Configuration file for M/EEG artefact detection
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_artefact.m 5592 2013-07-24 16:25:55Z vladimir $


%--------------------------------------------------------------------------
% D
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 1];
D.help   = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% mode
%--------------------------------------------------------------------------
mode = cfg_menu;
mode.tag = 'mode';
mode.name = 'Mode';
mode.labels = {'Reject', 'Mark'};
mode.val = {'reject'};
mode.values = {'reject', 'mark'};
mode.help = {'Action mode reject - to set trials and channels as bad',...
    'mark - just create artefact events, set channels as bad if mostly artefactual'};

%--------------------------------------------------------------------------
% badchanthresh
%--------------------------------------------------------------------------
badchanthresh         = cfg_entry;
badchanthresh.tag     = 'badchanthresh';
badchanthresh.name    = 'Bad channel threshold';
badchanthresh.strtype = 'r';
badchanthresh.num     = [1 1];
badchanthresh.val     = {0.2};
badchanthresh.help    = {'Fraction of trials with artefacts ', ...
    'above which an M/EEG channel is declared as bad.'};

%--------------------------------------------------------------------------
% methods
%--------------------------------------------------------------------------
artefact_funs = spm_select('List',spm('dir'),'^spm_eeg_artefact_.*\.m$');
artefact_funs = cellstr(artefact_funs);

fun      = cfg_choice;
fun.tag  = 'fun';
fun.name = 'Detection algorithm';
for i = 1:numel(artefact_funs)
    fun.values{i} = feval(spm_file(artefact_funs{i},'basename'));
end

methods      = cfg_branch;
methods.tag  = 'methods';
methods.name = 'Method';
methods.val  = {spm_cfg_eeg_channel_selector, fun};

%--------------------------------------------------------------------------
% methodsrep
%--------------------------------------------------------------------------
methodsrep        = cfg_repeat;
methodsrep.tag    = 'methodsrep';
methodsrep.name   = 'How to look for artefacts';
methodsrep.help   = {'Choose channels and methods for artefact detection'};
methodsrep.values = {methods};
methodsrep.num    = [1 Inf];

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''a''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'a'};

%--------------------------------------------------------------------------
% append
%--------------------------------------------------------------------------
append = cfg_menu;
append.tag = 'append';
append.name = 'Append';
append.labels = {'yes', 'no'};
append.val = {true};
append.values = {true, false};
append.help = {'Append new artefacts to already marked or overwrite'};

%--------------------------------------------------------------------------
% M/EEG Artefact detection
%--------------------------------------------------------------------------
artefact          = cfg_exbranch;
artefact.tag      = 'artefact';
artefact.name     = 'Artefact detection';
artefact.val      = {D, mode, badchanthresh, append, methodsrep, prefix};
artefact.help     = {'Detect artefacts in epoched M/EEG data.'};
artefact.prog     = @eeg_artefact;
artefact.vout     = @vout_eeg_artefact;
artefact.modality = {'EEG'};

%==========================================================================
% function out = eeg_artefact(job)
%==========================================================================
function out = eeg_artefact(job)
% construct the S struct
S.D = job.D{1};
S.mode = job.mode;
S.badchanthresh = job.badchanthresh;

for i = 1:numel(job.methods)
    S.methods(i).channels = spm_cfg_eeg_channel_selector(job.methods(i).channels);
    
    fun = fieldnames(job.methods(i).fun);
    fun = fun{1};
    
    S.methods(i).fun = fun;
    S.methods(i).settings = job.methods(i).fun.(fun);
end
    
out.D = spm_eeg_artefact(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

%==========================================================================
% function dep = vout_eeg_artefact(job)
%==========================================================================
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
