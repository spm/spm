D = spm_eeg_load('spm8_example.mat');

montage.labelorg = D.chanlabels;

montage.labelnew = [montage.labelorg(1:128), 'HEOG', 'VEOG'];

tra = eye(D.nchannels);
tra(131:end, :) = [];

% HEOG
tra(129, 129) = 0;
tra(129, [131 130]) = [1 -1];


% VEOG
tra(130, 130) = 0;
tra(130, [130 129]) = [1 -1];

montage.tra = tra;

save MONT_EXP montage
