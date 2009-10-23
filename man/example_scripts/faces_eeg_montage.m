D = spm_eeg_load('dspm8_faces_run1.mat');

montage.labelorg = D.chanlabels;

montage.labelnew = [montage.labelorg(1:128), 'HEOG', 'VEOG'];

tra = eye(D.nchannels);
tra(129:end, :) = [];
tra = detrend(tra, 'constant');

% HEOG
tra(129, [131 132]) = [1 -1];

% VEOG
tra(130, [135 136]) = [1 -1];

montage.tra = tra;

save faces_eeg_montage.mat montage
