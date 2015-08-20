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
% $Id: trialonset.m 1373 2008-04-11 14:24:03Z spm $

res = getset(this, 'trials', 'tag', varargin{:});
