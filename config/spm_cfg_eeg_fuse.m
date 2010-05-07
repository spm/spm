function S = spm_cfg_eeg_fuse
% configuration file for fusing M/EEG files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_fuse.m 3881 2010-05-07 21:02:57Z vladimir $

rev = '$Rev: 3881 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [2 Inf];
D.help = {'Select the M/EEG mat files.'};


S = cfg_exbranch;
S.tag = 'fuse';
S.name = 'M/EEG Fusion';
S.val = {D};
S.help = {'Fuse EEG/MEG data.'};
S.prog = @eeg_fuse;
S.vout = @vout_eeg_fuse;
S.modality = {'EEG'};

function out = eeg_fuse(job)
% construct the S struct
S.D = strvcat(job.D{:});

out.D = spm_eeg_fuse(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_fuse(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Fused Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Fused Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


