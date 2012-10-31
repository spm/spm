function res = reject(this, varargin)
% Method for getting/setting rejection flags
% FORMAT res = reject(this, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: reject.m 5025 2012-10-31 14:44:13Z vladimir $


warning_flexible('reject method is deprecated. Use ''badtrials'' instead');

res = getset(this, 'trials', 'bad', varargin{:});

if iscell(res)
    res = [res{:}];
end

