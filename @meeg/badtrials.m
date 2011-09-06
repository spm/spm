function res = badtrials(this, varargin)
% Method for getting/setting bad trials
% FORMAT res = badtrials(this)
% _______________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: badtrials.m 4465 2011-09-06 17:02:49Z guillaume $



if length(varargin) == 2 && ~isempty(varargin{1})
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 | (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

if numel(varargin) >= 1  && ~isempty(varargin{1})  
    if ~(all(varargin{1} >= 1) && all(varargin{1} <= ntrials(this)))
        error('Trial number out of range.');
    end
end

if numel(varargin) >= 2
    ubad = unique(varargin{2});
    if isempty(ubad) | ~all(ismember(ubad, [0 1]))
        error('Illegal bad flags (should be 0 or 1)');
    end
end

res = getset(this, 'trials', 'bad', varargin{:});


if isempty(varargin)
    if iscell(res)
        res = [res{:}];
    end
    
    res = find(res);
end