function res = badchannels(this, varargin)
% Method for getting/setting bad channels
% FORMAT res = badchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: badchannels.m 3942 2010-06-21 14:03:28Z vladimir $


if length(varargin) == 2
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 | (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

if numel(varargin) >= 1    
    if ~(varargin{1} >= 1 & varargin{1} <= nchannels(this))
        error('Channel number of out range.');
    end
end

if numel(varargin) >= 2
    ubad = unique(varargin{2});
    if isempty(ubad) | ~all(ismember(ubad, [0 1]))
        error('Illegal bad flags (should be 0 or 1)');
    end
end

res = getset(this, 'channels', 'bad', varargin{:});


if isempty(varargin)
    if iscell(res)
        res = [res{:}];
    end
    
    res = find(res);
end