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
% $Id: selectchannels.m 5675 2013-10-09 14:27:17Z vladimir $

if ischar(channels)
    channels = {channels};
end

chanind = [];

for i = 1:numel(channels)
    if ismember(upper(channels{i}), ...
            {'ALL', 'EOG', 'ECG', 'EMG', 'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR', 'MEGCOMB', 'REF', 'REFMAG', 'REFGRAD', 'LFP'})
        chanind = [chanind indchantype(this, upper(channels{i}))];
    elseif strncmpi('regexp_', channels{i}, 7)        
        re        = channels{i}(8:end);
        match     = regexp(chanlabels(this), re);
        chanind   = [chanind find(~cellfun('isempty', match))];
    else
        chanind = [chanind indchannel(this, channels{i})];
    end
    
    if any(size(chanind) == 0)
        chanind = [];
    end
end

chanind = unique(chanind);