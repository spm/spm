function res = time(this, ind, format)
% Method for getting the time axis
% FORMAT res = time(this, ind, format)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: time.m 5025 2012-10-31 14:44:13Z vladimir $

if this.Nsamples>0
    res = (0:(this.Nsamples-1))./this.Fsample + this.timeOnset;
else
    res = [];
end

if nargin>1 && ~isempty(ind)
    res = res(ind);
end

if nargin > 2
    if strcmp(format, 'ms')
        res = res*1000;
    end
end