function res = conditions(this, varargin)
% Method for getting condition labels, over trials
% FORMAT res = conditions(this, ind, conditionlabels)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: conditions.m 5025 2012-10-31 14:44:13Z vladimir $

res = getset(this, 'trials', 'label', varargin{:});

if nargin == 1 && ~iscell(res)
    res = {res};
end
