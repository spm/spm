function S = spm_cfg_eeg_merge
% configuration file for merging of M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_merge.m 1262 2008-03-28 09:28:42Z stefan $

D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [1 inf];
D.help = {'Select the M/EEG mat file.'};

S = cfg_exbranch;
S.tag = 'eeg_merge';
S.name = 'M/EEG merging';
S.val = {D};
S.help = {'Merge EEG/MEG data.'};
S.prog = @eeg_merge;
S.vout = @vout_eeg_merge;
S.modality = {'EEG'};


function out = eeg_merge(job)
% construct the S struct
S.D = strvcat(job.D);
out.D = spm_eeg_merge(S);

function dep = vout_eeg_merge(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Merge Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


