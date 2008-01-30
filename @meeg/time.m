function res = time(obj)
% Method for getting the time axis
% FORMAT res = time(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: time.m 1125 2008-01-30 12:12:18Z vladimir $

if obj.Nsamples>0
    res = (0:(obj.Nsamples-1))./obj.Fsample + obj.timeOnset;
else
    res = [];
end