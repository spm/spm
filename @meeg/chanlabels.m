function res = chanlabels(obj, ind)
% Method for getting the channel labels
% FORMAT res = chanlabels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chanlabels.m 1239 2008-03-25 15:28:16Z vladimir $

if obj.Nsamples>0
    res = {obj.channels.label};
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

% res = strvcat(res);
