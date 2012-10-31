function res = events(this, varargin)
% Method for getting/setting events per trial
% FORMAT res = events(this, ind, event)
%   ind = indices of trials
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: events.m 5025 2012-10-31 14:44:13Z vladimir $

res = getset(this, 'trials', 'events', varargin{:});
