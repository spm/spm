function res = badchannels(this, varargin)
% Method for getting/setting bad channels
% FORMAT res = badchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: badchannels.m 5076 2012-11-23 16:05:21Z vladimir $

if length(varargin) == 2 && ~isempty(varargin{1})
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 || (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

if numel(varargin) >= 1  && ~isempty(varargin{1})  
    if ~(all(varargin{1} >= 1) && all(varargin{1} <= nchannels(this)))
        error('Channel number out of range.');
    end
end

if numel(varargin) >= 2 && ~isempty(varargin{1})
    ubad = unique(varargin{2});
    if isempty(ubad) || ~all(ismember(ubad, [0 1]))
        error('Illegal bad flags (should be 0 or 1)');
    end
end

if this.montage.Mind == 0
    res = getset(this, 'channels', 'bad', varargin{:});
else 
    res = getset(this.montage.M(this.montage.Mind), 'channels', 'bad', varargin{:});
end

% Return channel indices if called without arguments and [0, 1] if called
if numel(varargin) <= 1 % get
    if iscell(res)
        res = [res{:}];
    end
    if isempty(varargin)
        res = find(res);
    end
end
    
