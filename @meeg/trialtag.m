function res = trialtag(this, varargin)
% Method for getting/setting trial tag
% FORMAT res = trialtag(this, ind, tag)
%   ind = indices of trials
% The user can put any data here that will be attached to
% the respective trial. This is useful e.g. to make sure the
% relation between regressors and data is not broken when
% removing bad trials or merging files.
% _______________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: trialtag.m 6618 2015-12-01 16:25:38Z spm $

res = getset(this, 'trials', 'tag', varargin{:});
