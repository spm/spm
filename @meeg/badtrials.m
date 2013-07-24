function res = badtrials(this, varargin)
% Method for getting/setting bad trials
% FORMAT res = badtrials(this)
% _______________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: badtrials.m 5592 2013-07-24 16:25:55Z vladimir $


if length(varargin) == 2 && ~isempty(varargin{1})
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 || (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

if numel(varargin) >= 2
    ubad = unique(varargin{2});
    if isempty(ubad) || ~all(ismember(ubad, [0 1]))
        error('Illegal bad flags (should be 0 or 1)');
    end
end

res = getset(this, 'trials', 'bad', varargin{:});


% Return trial indices if called without arguments and [0, 1] if called
if numel(varargin) <= 1 % get
    if iscell(res)
        res = [res{:}];
    end
    res = logical(res);
    if isempty(varargin)
        res = find(res);
    end
end
