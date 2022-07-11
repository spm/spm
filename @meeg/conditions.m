function res = conditions(this, varargin)
% Method for getting condition labels, over trials
% FORMAT res = conditions(this, ind, conditionlabels)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = getset(this, 'trials', 'label', varargin{:});

if nargin == 1 && ~iscell(res)
    res = {res};
end
