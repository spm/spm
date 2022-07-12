% These are some examples of how to convert data between SPM and FieldTrip
% Tested on SPM EEG MMN example with headmodel specified.

% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging

spm('defaults', 'eeg');

D = spm_eeg_load('maeafdfMspmeeg_subject1.mat');

% Exporting a data subset to FT raw struct
raw = D.ftraw(D.indchantype('EEG'), D.indsample(-0.05):D.indsample(0.1), D.indtrial('std'));

% Exporting a data subset to FT timelock struct
timelock = D.fttimelock(D.indchantype('EEG'), D.indsample(-0.05):D.indsample(0.1), D.indtrial('std'));

% Importing from SPM to FT
hdr = ft_read_header('maeafdfMspmeeg_subject1.mat');

cfg = [];
cfg.dataset = 'maeafdfMspmeeg_subject1.mat';
data = ft_preprocessing(cfg);


% Importing from FT to SPM
Dnew = spm_eeg_ft2spm(data, 'ft2spm_example');

% Exporting the head model for use in FT.
% This might not work with data moved to a different computer
[vol_sens] = spm_eeg_inv_get_vol_sens(D);
headmodel = vol_sens.EEG.vol;
sens =  vol_sens.EEG.sens;
