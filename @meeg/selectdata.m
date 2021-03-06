function res = selectdata(this, varargin)
% Selects data using channel labels, time and condition labels as indices
% FORMAT res = selectdata(D, chanlabel, timeborders, condition)
%        res = selectdata(D, chanlabel, freqborders, timeborders, condition)
%
%  D - meeg object
%  chanlabel   - channel label, cell array of labels or [] (for all channels)
%  timeborders - [start end] in sec or [] for all times
%  freqborders - [start end] in Hz or [] for all frequencis (for TF datasets only)
%  condition   - condition label, cell array of labels or [] (for all conditions)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if this.Nsamples == 0
    res = [];
    return;
end

if (isequal(transformtype(this), 'time') && numel(varargin)<3) ||...
        (strncmpi(transformtype(this),'TF',2) &&  numel(varargin)<4)
    error('Insufficient number of arguments for data selection');
end

chanlabel = varargin{1};

if isequal(transformtype(this), 'time')
    timeborders = varargin{2};
    condition   = varargin{3};
elseif strncmpi(transformtype(this),'TF',2)
    freqborders = varargin{2};
    timeborders = varargin{3};
    condition   = varargin{4};

    if isempty(freqborders)
        freqind = 1:nfrequencies(this);
    else
        freqind = indfrequency(this, freqborders(1)):indfrequency(this, freqborders(2));
    end
else
    error('Unsupported transform type.');
end

if isempty(chanlabel)
    chanind = 1:nchannels(this);
else
    chanind = indchannel(this, chanlabel);
end

if isempty(timeborders)
    timeind = 1:nsamples(this);
else
    timeind = indsample(this, timeborders(1)):indsample(this, timeborders(2));
end

if isempty(condition)
    trialind = 1:ntrials(this);
else
    trialind = indtrial(this, condition);
end

ss.type = '()';
if isequal(transformtype(this), 'time')
    ss.subs = {chanind, timeind, trialind};
    res = subsref(this,ss);
else
    ss.subs = {chanind, freqind, timeind, trialind};
    res = subsref(this,ss);
end
