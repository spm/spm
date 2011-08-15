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
% $Id: chantype.m 4432 2011-08-15 12:43:44Z christophe $

if this.montage.Mind==0
    res = getset(this, 'channels', 'type', varargin{:});
else
    res = getset(this.montage.M(this.montage.Mind), 'channels', 'type', varargin{:});
end
