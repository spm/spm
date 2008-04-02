function S = spm_cfg_eeg_grandmean
% configuration file for averaging evoked responses
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_grandmean.m 1292 2008-04-02 14:17:31Z volkmar $

rev = '$Rev';
D = cfg_files;
D.tag = 'D';
D.name = 'File Names';
D.filter = 'mat';
D.num = [1 inf];
D.help = {'Select the M/EEG mat file.'};

Pout = cfg_entry;
Pout.tag = 'Pout';
Pout.name = 'Output filename';
Pout.strtype = 's';
Pout.num = [1 inf];
Pout.help = {'Choose filename'};

S = cfg_exbranch;
S.tag = 'eeg_grandmean';
S.name = 'M/EEG grandmean';
S.val = {D Pout};
S.help = {'Average multiple evoked responses'};
S.prog = @eeg_grandmean;
S.vout = @vout_eeg_grandmean;
S.modality = {'EEG'};


function out = eeg_grandmean(job)
% construct the S struct
S.P = strvcat(job.D);
S.Pout = job.Pout;

out.D = spm_eeg_grandmean(S);

function dep = vout_eeg_grandmean(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Take grandmean of evoked responses';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


