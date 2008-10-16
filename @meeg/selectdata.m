function res = selectdata(this, chanlabel, timeborders, condition)
% Selects data using channel labels, time and condition labels as indices
% FORMAT res = selectdata(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectdata.m 2347 2008-10-16 12:58:33Z vladimir $

if this.Nsamples == 0
    res = [];
    return;
end

if nargin<4
    error('Insufficient number of arguments for data selection');
end

if isempty(chanlabel)
    chanind = 1:nchannels(this);
else
    if ischar(chanlabel)
        chanlabel = {chanlabel};
    end
    [junk, chanind] = spm_match_str(chanlabel, chanlabels(this));
end

if isempty(timeborders)
    timeind = 1:nsamples(this);
else
    timeAxis = time(this);
    timeStart = min(find(timeAxis>=timeborders(1)));
    timeEnd = max(find(timeAxis<=timeborders(2)));
    timeind = timeStart:timeEnd;
end


if isempty(condition)
    trialind = 1:ntrials(this);
else
    if ischar(condition)
        condition = {condition};
    end
    [junk, trialind] = spm_match_str(condition, {this.trials.label});
end

res = double(this.data.y(chanind, timeind, trialind));