function res = chanlabels(obj, ind)
% Method for getting the channel labels
% FORMAT res = chanlabels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chanlabels.m 1236 2008-03-20 18:15:33Z stefan $

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

res = strvcat(res);
