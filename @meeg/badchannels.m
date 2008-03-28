function res = badchannels(this, varargin)
% Method for getting/setting bad channels 
% FORMAT res = badchannels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$

    
res = getset(this, 'channels', 'bad', varargin{:});

if isempty(varargin)
    res = find(res);
end