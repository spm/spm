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
% $Id: chantype.m 5933 2014-03-28 13:22:28Z vladimir $

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
