function res = conditions(this, varargin)
% Method for getting condition labels, over trials
% FORMAT res = conditions(this, ind, conditionlabels)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
%$Id $

res = getset(this, 'trials', 'label', varargin{:});