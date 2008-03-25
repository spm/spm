function res = badchannels(this, ind, flag)
% Method for getting/setting bad channels 
% FORMAT res = badchannels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$

switch nargin
    case 1
        res = getbadchannels(this);        
    case 3
        res = setbadchannels(this, ind, flag);
    otherwise
end

function res = getbadchannels(this)
res = find(cat(1, this.channels(:).bad));

function this = setbadchannels(this, ind, flag)
[this.channels(ind).bad] = deal(flag);