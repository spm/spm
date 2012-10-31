function res = indtrial(this, label)
% Method for getting trial indices based on condition labels
% FORMAT res = indtrial(this, label)
% this       - MEEG object
% label      - string or cell array of labels, 'GOOD' and 'BAD'
%              can be added to list of labels to select only
%              good or bad trials respectively
% res        - vector of trial indices matching condition labels
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indtrial.m 5025 2012-10-31 14:44:13Z vladimir $

if ischar(label)
    label = {label};
end

[junk, res] = match_str(label, conditions(this));

if ismember('GOOD', upper(label)) && ~ismember('BAD', upper(label))
    res = setdiff(res, badtrials(this));
elseif ismember('BAD', upper(label)) && ~ismember('GOOD', upper(label))
    res = intersect(res, badtrials(this));
end