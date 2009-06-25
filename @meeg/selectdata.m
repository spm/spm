function res = selectdata(this, chanlabel, timeborders, condition)
% Selects data using channel labels, time and condition labels as indices
% FORMAT res = selectdata(D, chanlabel, timeborders, condition)
%
%  D - meeg object
%  chanlabel - channel label, cell array of labels or [] (for all channels)
%  timeborders - [start end] in sec or [] for all times
%  condition   - condition label, cell array of labels or [] (for all conditions)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectdata.m 3219 2009-06-25 10:15:36Z vladimir $

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