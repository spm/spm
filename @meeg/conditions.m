function res = conditions(this, varargin)
% Method for getting condition labels, over trials
% FORMAT res = conditions(this, ind, conditionlabels)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: conditions.m 1373 2008-04-11 14:24:03Z spm $

res = getset(this, 'trials', 'label', varargin{:});

if nargin ==1 & ~iscell(res)
    res = {res};
end
