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
% $Id: chantype.m 2105 2008-09-17 15:34:09Z vladimir $


if length(varargin)>=2
    ind = varargin{1};
    type = varargin{2};

    if isempty(ind)
        ind = 1:nchannels(this);
    end

    if isempty(type)
        types = {'EEG', 'MEG', 'MEGREF', 'EMG', 'EOG'};

        for i=1:length(types)
            if isempty(ind) break; end
            foundind = spm_match_str(chanlabels(this, ind), ...
                ft_channelselection(types(i), chanlabels(this, ind)));
            if ~isempty(foundind)
                this = chantype(this, ind(foundind), types{i});
                ind = setdiff(ind, ind(foundind));
            end
        end

        if isfield(this, 'origchantypes')
            ind1 = setdiff(ind, strmatch('unknown', this.other.origchantypes, 'exact'));
            if ~isempty(ind1)
                this = chantype(this, ind1, upper(this.other.origchantypes(ind1)));
                ind = setdiff(ind, ind1);
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