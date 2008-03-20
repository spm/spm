function res = time(obj, ind, format)
% Method for getting the time axis
% FORMAT res = time(obj, ind, format)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: time.m 1236 2008-03-20 18:15:33Z stefan $

if obj.Nsamples>0
    res = (0:(obj.Nsamples-1))./obj.Fsample + obj.timeOnset;
else
    res = [];
end

if exist('ind') == 1
    res = res(ind);
end

if exist('format') ==1
    if strcmp(format, 'ms')
        res = res*1000;
    end
end