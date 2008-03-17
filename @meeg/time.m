function res = time(this)
% Method for getting the time axis
% FORMAT res = time(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: time.m 1219 2008-03-17 17:35:12Z vladimir $

if this.Nsamples>0
    res = (0:(this.Nsamples-1))./this.Fsample + this.timeOnset;
else
    res = [];
end