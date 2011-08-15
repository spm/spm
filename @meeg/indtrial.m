function res = indtrial(this, label)
% Method for getting trial indices based on channel labels
% FORMAT  res = indtrial(this, label)
% this       - MEEG object
% label      - string or cell array of labels
%
% res        - vector of trial indices matching condition labels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indtrial.m 4432 2011-08-15 12:43:44Z christophe $

if ischar(label)
    label = {label};
end

[junk, res] = match_str(label, conditions(this));