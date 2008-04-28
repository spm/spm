function res = reject(this, varargin)
% Method for getting/setting rejection flags
% FORMAT res = reject(this, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: reject.m 1490 2008-04-28 11:16:29Z vladimir $


res = getset(this, 'trials', 'bad', varargin{:});

