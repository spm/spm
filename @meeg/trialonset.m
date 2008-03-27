function res = trialonset(this, varargin)
% Method for getting/setting trial onset times
% FORMAT res = trialonset(this, ind, onset)
%   ind = indices of trials
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $

res = getset(this, 'trials', 'onset', varargin{:});