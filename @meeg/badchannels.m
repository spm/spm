function res = badchannels(this, varargin)
% Method for getting/setting bad channels 
% FORMAT res = badchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: badchannels.m 1373 2008-04-11 14:24:03Z spm $

    
if length(varargin) == 2
    % make sure that the two inputs for set are the same length
    if ~(length(varargin{2}) == 1 | (length(varargin{1}) == length(varargin{2})))
        error('Use either same vector length or scalar for value');
    end
end

res = getset(this, 'channels', 'bad', varargin{:});

if isempty(varargin)
    res = find(res);
end