function res = indchannel(this, label)
% Method for getting channel indices based on channel labels
% FORMAT  res = indchannel(this, label)
% this       - MEEG object
% label      - string or cell array of labels
%
% res        - vector of channel indices matching labels
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ischar(label)
    label = {label};
end

[junk, res] = match_str(label, chanlabels(this));
