function S = spm_cfg_eeg_average
% configuration file for M/EEG epoching
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_average.m 2870 2009-03-12 13:43:35Z guillaume $

rev = '$Rev: 2870 $';
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

S = cfg_exbranch;
S.tag = 'eeg_average';
S.name = 'M/EEG averaging';
S.val = {D};
S.help = {'Average epoched EEG/MEG data.'};
S.prog = @eeg_average;
S.vout = @vout_eeg_average;
S.modality = {'EEG'};


function out = eeg_average(job)
% construct the S struct
S.D = job.D{1};
S.review = false;

out.D = spm_eeg_average(S);
out.Dfname = {out.D.fname};

function dep = vout_eeg_average(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Average Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Averaged Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


