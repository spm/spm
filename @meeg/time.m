function res = time(this, ind, format)
% Method for getting the time axis
% FORMAT res = time(this, ind, format)
%__________________________________________________________________________

% Vladimir Litvak, Stefan Kiebel
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


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
