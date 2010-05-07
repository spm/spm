function S = spm_cfg_eeg_merge
% configuration file for merging of M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Volkmar Glauche
% $Id: spm_cfg_eeg_merge.m 3881 2010-05-07 21:02:57Z vladimir $

rev = '$Rev: 3881 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [2 Inf];
D.help = {'Select the M/EEG mat files.'};

file = cfg_entry;
file.tag = 'file';
file.name = 'Files to which the rule applies';
file.strtype = 's';
file.val = {'.*'};
file.help = {'Regular expression to match the files to which the rule applies (default - all)'}; 

labelorg = cfg_entry;
labelorg.tag = 'labelorg';
labelorg.name = 'Original labels to which the rule applies';
labelorg.strtype = 's';
labelorg.val = {'.*'};
labelorg.help = {'Regular expression to match the original condition labels to which the rule applies (default - all)'}; 

labelnew = cfg_entry;
labelnew.tag = 'labelnew';
labelnew.name = 'New label for the merged file';
labelnew.strtype = 's';
labelnew.val = {'#labelorg#'};
labelnew.help = {['New condition label for the merged file. Special tokens can be used as part of the name. '...
    '#file# will be replaced by the name of the original file, #labelorg# will be replaced by the original '...
    'condition labels.']}; 

rule = cfg_branch;
rule.tag = 'rule';
rule.name = 'Recoding rule';
rule.val  = {file, labelorg, labelnew};
rule.help = {'Recoding rule. The default means that all trials will keep their original label.'}; 


rules = cfg_repeat;
rules.tag = 'unused';
rules.name = 'Condition label recoding rules';
rules.values = {rule};
rules.num = [1 Inf];
rules.val = {rule};
rules.help = {['Specify the rules for translating condition labels from ' ...
    'the original files to the merged file. Multiple rules can be specified. The later ' ...
    'rules have precedence. Trials not matched by any of the rules will keep their original labels.']};

S = cfg_exbranch;
S.tag = 'merge';
S.name = 'M/EEG Merging';
S.val = {D, rules};
S.help = {'Merge EEG/MEG data.'};
S.prog = @eeg_merge;
S.vout = @vout_eeg_merge;
S.modality = {'EEG'};

function out = eeg_merge(job)
% construct the S struct
S.D = strvcat(job.D{:});
S.recode = job.rule;

out.D = spm_eeg_merge(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_merge(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Merged Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Merged Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


