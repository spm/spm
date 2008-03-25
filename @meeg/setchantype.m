function this = setchantype(this, ind, type)
% Method for getting the time axis
% FORMAT setchantype(this, ind, type)
%   ind - channel index
%   type - type (string: 'EEG', 'MEG', 'LFP' etc.)
%
% FORMAT setchantype(this, ind), setchantype(this)
% Sets channel types to default using Fieldtrip channelselection
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: setchantype.m 1240 2008-03-25 16:34:30Z vladimir $

if nargin == 1 || isempty(ind)
    ind = 1:nchannels(this);
end

if nargin == 3
    [this.channels(ind).type] = deal(type);
    return
end

if exist('channelselection') == 2
    eeg_types = {'EEG', 'EEG1020', 'EEG1010', 'EEG1005', 'EEGBHAM', 'EEGREF'};
    other_types = {'MEG', 'EMG', 'EOG'};

    for i=1:length(eeg_types)
        foundind = spm_match_str(chanlabels(this, ind), channelselection(eeg_types(i), chanlabels(this, ind)));
        if ~isempty(foundind)
            this = setchantype(this, ind(foundind), 'EEG');
            ind = setdiff(ind, ind(foundind));
        end
    end

    for i=1:length(other_types)
        foundind = spm_match_str(chanlabels(this, ind), channelselection(other_types(i), chanlabels(this, ind)));
        if ~isempty(foundind)
            this = setchantype(this, ind(foundind), other_types{i});
            ind = setdiff(ind, ind(foundind));
        end
    end
else
    warning('Fieldtrip not available, setting all types to ''other''');
end

if ~isempty(ind)
    this = setchantype(this, ind, 'Other');
end
