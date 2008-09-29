function S = spm_cfg_eeg_merge
% configuration file for merging of M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Volkmar Glauche
% $Id: spm_cfg_eeg_merge.m 2225 2008-09-29 12:25:27Z stefan $

rev = '$Rev: 2225 $';
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
files.num = [2 inf];
files.values = {file};

S = cfg_exbranch;
S.tag = 'eeg_merge';
S.name = 'M/EEG Merging';
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
out.Dfname = {out.D.fname};

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


