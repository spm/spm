function res = chanlabels(obj)
% Method for getting the channel labels
% FORMAT res = chanlabels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chanlabels.m 1125 2008-01-30 12:12:18Z vladimir $

if obj.Nsamples>0
    res = {obj.channels.label};
else
    res = [];
end