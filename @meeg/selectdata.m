function res = selectdata(obj, chanlabel, timeborders, condition)
% Selects data using channel labels, time and condition labels as indices
% FORMAT res = selectdata(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectdata.m 1239 2008-03-25 15:28:16Z vladimir $

if obj.Nsamples == 0
    res = [];
    return;
end

if nargin<4
    error('Insufficient number of arguments for data selection');
end

if isempty(chanlabel)
    chanind = 1:nchannels(obj);
else
    if ischar(chanlabel)
        chanlabel = {chanlabel};
    end
    [junk, chanind] = spm_match_str(chanlabel, chanlabels(obj));
end

if isempty(timeborders)
    timeind = 1:nsamples(obj);
else
    timeAxis = time(obj);
    timeStart = min(find(timeAxis>=timeborders(1)));
    timeEnd = max(find(timeAxis<=timeborders(2)));
    timeind = timeStart:timeEnd;
end


if isempty(condition)
    trialind = 1:ntrials(obj);
else
    if ischar(condition)
        condition = {condition};
    end
    [junk, trialind] = spm_match_str(condition, {obj.trials.label});
end

res = obj.data.y(chanind, timeind, trialind);