function S = spm_cfg_eeg_merge
% configuration file for merging of M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Volkmar Glauche
% $Id: spm_cfg_eeg_merge.m 1648 2008-05-15 09:49:36Z stefan $

rev = '$Rev: 1648 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

string = cfg_entry;
string.tag = 'str';
string.name = 'Label';
string.strtype = 's';
string.num = [1 Inf];

strings = cfg_repeat;
strings.tag = 'unused';
strings.name = 'Labels';
strings.values = {string};

file = cfg_branch;
file.tag = 'file';
file.name = 'File info';
file.val  = {D strings};

files = cfg_repeat;
files.tag = 'unused';
files.name = 'Files';
files.values = {file};

S = cfg_exbranch;
S.tag = 'eeg_merge';
S.name = 'M/EEG merging';
S.val = {files};
S.help = {'Merge EEG/MEG data.'};
S.prog = @eeg_merge;
S.vout = @vout_eeg_merge;
S.modality = {'EEG'};


function out = eeg_merge(job)
% construct the S struct
S.D = strvcat(cat(1,job.file(:).D));
for i = 1:length(job.file)
    for j = 1:numel(job.file(i).str)
        S.recode{i}{j} = job.file(i).str{j};
    end
end

out.D = spm_eeg_merge(S);

function dep = vout_eeg_merge(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Merge Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


