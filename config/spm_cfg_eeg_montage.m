function S = spm_cfg_eeg_montage
% configuration file for reading montage files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_montage.m 4212 2011-02-23 17:50:55Z vladimir $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

montage = cfg_files;
montage.tag = 'montage';
montage.name = 'Montage file name';
montage.filter = 'mat';
montage.num = [1 1];
montage.help = {'Select a montage file.'};

keepothers = cfg_menu;
keepothers.tag = 'keepothers';
keepothers.name = 'Keep other channels';
keepothers.labels = {'Yes', 'No'};
keepothers.values = {'yes', 'no'};
keepothers.val = {'no'};
keepothers.help = {'Specify whether you want to keep channels that are not contributing to the new channels'};

S = cfg_exbranch;
S.tag = 'montage';
S.name = 'M/EEG Montage';
S.val = {D montage keepothers};
S.help = {'Apply a montage (linear transformation) to EEG/MEG data.'};
S.prog = @eeg_montage;
S.vout = @vout_eeg_montage;
S.modality = {'EEG'};

function out = eeg_montage(job)
% construct the S struct
S.D = job.D{1};
S.montage = job.montage{1};
S.keepothers = job.keepothers;
out.D = spm_eeg_montage(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_montage(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'montaged data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Montaged Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


