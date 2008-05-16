function res = badchannels(this, varargin)
% Method for getting/setting bad channels 
% FORMAT res = badchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: badchannels.m 1673 2008-05-16 15:32:22Z vladimir $

    
if length(varargin) == 2
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 | (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

if length(varargin) >= 1
    if ~(varargin{1} >= 1 & varargin{1} <= nchannels(this))
        error('Channel number of out range.');
    end
end

res = getset(this, 'channels', 'bad', varargin{:});


if isempty(varargin)
    if iscell(res)
      res = [res{:}];
    end
    
    res = find(res);
end