function chanind = selectchannels(this, channels)
% Method for getting channel indices based on labels and/or types
% FORMAT  res = selectchannels(this, label)
% this       - MEEG object
% channels   - string or cell array of labels that may also include 
%              'all', or types ('EEG', 'MEG' etc.)
%
% res        - vector of channel indices matching labels
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectchannels.m 5212 2013-01-26 13:16:36Z vladimir $

if ischar(channels)
    channels = {channels};
end

chanind = [];

for i = 1:numel(channels)
    if ismember(upper(channels{i}), ...
            {'ALL', 'EOG', 'ECG', 'EMG', 'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR', 'MEGCOMB', 'REF', 'REFMAG', 'REFGRAD', 'LFP'})
        chanind = [chanind indchantype(this, upper(channels{i}))];
    else
        chanind = [chanind indchannel(this, channels{i})];
    end
end

chanind = unique(chanind);