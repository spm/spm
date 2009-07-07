function res = indtrial(this, label)
% Method for getting channel indices based on channel labels
% FORMAT  res = indtrial(this, label)
% this       - MEEG object
% label      - string or cell array of labels
%
% res        - vector of trial indices matching condition labels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indtrial.m 3254 2009-07-07 15:18:54Z vladimir $

if ischar(label)
    label = {label};
end

[junk, res] = match_str(label, conditions(this));