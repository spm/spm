function dictionary = spm_eeg_spmft_chan_dictionary
% Return a table of corresponce between SPM and FieldTrip channel types
% FORMAT dictionary = spm_eeg_spmft_chan_dictionary
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging

dictionary = {
            'eog',           'EOG';
            'eeg',           'EEG';
            'ecg',           'ECG';
            'lfp',           'LFP';
            'emg',           'EMG';
            'meg',           'MEG';
            'ref',           'REF';
            'megref'         'REF';
            'megmag',        'MEGMAG';
            'megplanar',     'MEGPLANAR';
            'meggrad',       'MEGGRAD';
            'refmag',        'REFMAG';
            'refgrad',       'REFGRAD';
            'refplanar'      'REFPLANAR';
            'unknown'        'Other';
            };
