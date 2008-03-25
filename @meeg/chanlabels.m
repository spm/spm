function res = chanlabels(this, ind)
% Method for getting/setting the channel labels
% FORMAT res = chanlabels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

switch nargin
    case 1
        res = getlabels(this);        
    case 2
        res = getlabels(this, ind);
    case 3
        res = setlabels(this, ind, name);
    otherwise
end

function res = getlabels(this, ind)
if this.Nsamples>0
    res = {this.channels.label};
else
    res = [];
end

if nargin > 1
    if all(ind > 0 & ind <= nchannels(obj))
        res = res(ind);
    else
        error('Indexed channels do not exist.');
    end
end


function this = setlabels(this, ind, name)
if iscell(label)
    [this.trials(ind).label] = deal(name);
else
    [this.trials(ind).label] = deal(cellstr(name));
end