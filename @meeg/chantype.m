function res = chantype(this, varargin)
% Method for setting/getting channel types
% FORMAT chantype(this, ind, type)
%   ind - channel index
%   type - type (string: 'EEG', 'MEG', 'LFP' etc.)
%
% FORMAT chantype(this, ind), chantype(this)
% Sets channel types to default using Fieldtrip channelselection
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chantype.m 1334 2008-04-09 14:28:11Z vladimir $


if length(varargin)>=2
    ind = varargin{1};
    type = varargin{2};

    if isempty(ind)
        ind = 1:nchannels(this);
    end

    if isempty(type)
        eeg_types = {'EEG', 'EEG1020', 'EEG1010', 'EEG1005', 'EEGBHAM', 'EEGREF'};
        other_types = {'MEG', 'EMG', 'EOG'};

        for i=1:length(eeg_types)
            if isempty(ind) break; end
            foundind = spm_match_str(chanlabels(this, ind), ft_channelselection(eeg_types(i), chanlabels(this, ind)));
            if ~isempty(foundind)
                this = chantype(this, ind(foundind), 'EEG');
                ind = setdiff(ind, ind(foundind));
            end
        end

        for i=1:length(other_types)
            if isempty(ind) break; end
            foundind = spm_match_str(chanlabels(this, ind), ft_channelselection(other_types(i), chanlabels(this, ind)));
            if ~isempty(foundind)
                this = chantype(this, ind(foundind), other_types{i});
                ind = setdiff(ind, ind(foundind));
            end
        end
        
        if ~isempty(ind)
            this = chantype(this, ind, 'Other');
        end

        res = this;
        return
    end
end

res = getset(this, 'channels', 'type', varargin{:});