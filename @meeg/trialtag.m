function res = trialtag(this, varargin)
% Method for getting/setting trial tag
% FORMAT res = trialtag(this, ind, tag)
%   ind = indices of trials
% The user can put any data here that will be attached to
% the respective trial. This is useful e.g. to make sure the
% relation between regressors and data is not broken when
% removing bad trials or merging files.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = getset(this, 'trials', 'tag', varargin{:});
