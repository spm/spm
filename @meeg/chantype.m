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
% $Id: chantype.m 2404 2008-10-27 17:34:07Z vladimir $


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

        if ~isempty(ind) && isfield(this, 'origchantypes')
            ind1 = setdiff(1:numel(this.other.origchantypes.label),...
                strmatch('unknown', this.other.origchantypes.type, 'exact'));
            
            [sel1, sel2] = spm_match_str(chanlabels(this, ind), this.other.origchantypes.label(ind1));
            
            if ~isempty(sel1)
                this = chantype(this, ind(sel1), upper(this.other.origchantypes.type(ind1(sel2))));
                ind(sel1) = [];
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