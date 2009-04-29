function S = spm_cfg_eeg_tf_rescale
% configuration file for rescaling spectrograms
%_______________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_cfg_eeg_tf_rescale.m 3090 2009-04-29 17:30:05Z will $

Fname = cfg_files;
Fname.tag = 'Fname';
Fname.name = 'File Names';
Fname.filter = 'mat';
Fname.num = [1 inf];
Fname.help = {'Select the M/EEG mat file.'};


S = cfg_exbranch;
S.tag = 'eeg_tf_rescale';
S.name = 'M/EEG TF Rescale';
S.val = {Fname};
S.help = {'Rescale spectrogram using eg Log Ratio operator'};
S.prog = @eeg_tf_rescale;
S.modality = {'EEG'};


function out = eeg_tf_rescale(job)
% construct the S struct
S = job;
S.Fname = strvcat(job.Fname);

