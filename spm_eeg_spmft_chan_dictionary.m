function dictionary = spm_eeg_spmft_chan_dictionary
% Returns a table of corresponce between SPM and FieldTrip channel types
%______________________________________________________________________
% Copyright (C) 2021 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spmft_chan_dictionary.m 8130 2021-08-05 13:15:12Z vladimir $

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