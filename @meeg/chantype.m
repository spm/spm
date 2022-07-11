function res = chantype(this, varargin)
% Method for setting/getting channel types
% FORMAT chantype(this, ind, type)
%   ind - channel index
%   type - type (string: 'EEG', 'MEG', 'LFP' etc.)
%
% FORMAT chantype(this, ind), chantype(this)
% Sets channel types to default using Fieldtrip channelselection
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if this.montage.Mind==0
    res = getset(this, 'channels', 'type', varargin{:});
else
    if nargin == 3
        this.montage.M(this.montage.Mind) = getset(this.montage.M(this.montage.Mind), 'channels', 'type', varargin{:});
        res = this;
    else
        res = getset(this.montage.M(this.montage.Mind), 'channels', 'type', varargin{:});
    end
end
