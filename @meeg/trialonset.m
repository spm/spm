function res = trialonset(this, varargin)
% Method for getting/setting trial onset times
% FORMAT res = trialonset(this, ind, onset)
%   ind = indices of trials
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = getset(this, 'trials', 'onset', varargin{:});
