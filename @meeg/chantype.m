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
% $Id: chantype.m 2078 2008-09-10 17:33:21Z vladimir $


if length(varargin)>=2
    ind = varargin{1};
    type = varargin{2};

    if isempty(ind)
        ind = 1:nchannels(this);
    end

    if isempty(type)
        eeg_types = {'EEG'};
        other_types = {'MEG', 'EMG', 'EOG'};

        for i=1:length(eeg_types)
            if isempty(ind) break; end
            if length(ind) == 1
                foundind = spm_match_str({chanlabels(this, ind)}, ...
                    ft_channelselection(eeg_types(i), {chanlabels(this, ind)}));
            else
                foundind = spm_match_str(chanlabels(this, ind), ...
                    ft_channelselection(eeg_types(i), chanlabels(this, ind)));
            end
            if ~isempty(foundind)
                this = chantype(this, ind(foundind), 'EEG');
                ind = setdiff(ind, ind(foundind));
            end
        end

        for i=1:length(other_types)
            if isempty(ind) break; end
            if length(ind) == 1
                foundind = spm_match_str({chanlabels(this, ind)}, ...
                    ft_channelselection(other_types(i), {chanlabels(this, ind)}));
            else
                foundind = spm_match_str(chanlabels(this, ind), ...
                    ft_channelselection(other_types(i), chanlabels(this, ind)));
            end
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