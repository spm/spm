function timelock = fttimelock(this, chanind, timeind, trialind, freqind)
% Method for converting meeg object to Fieldtrip timelock/freq struct
% FORMAT timelock = fttimelock(this, chanind, timeind, trialind, freqind)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fttimelock.m 5025 2012-10-31 14:44:13Z vladimir $

if ~islinked(this)
    error('There is no data linked to the object');
end

if nargin < 2 || isempty(chanind)
    chanind = 1:nchannels(this);
end

if nargin < 3 || isempty(timeind)
    timeind = 1:nsamples(this);
end

if nargin < 4 || isempty(trialind)
    trialind = 1:ntrials(this);
end

if strncmpi(transformtype(this),'TF',2) && ...
        (nargin < 5 || isempty(freqind))
     freqind = 1:nfrequencies(this);
end

timelock             = [];
timelock.label       = chanlabels(this, chanind)';

if isequal(transformtype(this), 'time')
    timelock.dimord  = 'rpt_chan_time';
    timelock.trial   =  permute(this.data(chanind, timeind, trialind), [3 1 2]);
    timelock.dimord  = 'rpt_chan_time';
elseif strncmpi(transformtype(this),'TF',2)
    timelock.dimord  = 'rpt_chan_freq_time';
    timelock.powspctrm = permute(this.data(chanind, freqind, timeind, trialind), [4 1 2 3]);
    timelock.freq    =  frequencies(this, freqind);
else
  error('Unknown transform type.');
end

timelock.time       = time(this, timeind);

clist      =  condlist(this);

condlabels = conditions(this, trialind);

timelock.trialinfo = 0*trialind;
for k = 1:numel(clist)
  fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
  sel = strcmp(clist{k}, condlabels);
  timelock.trialinfo(sel) = k;
end

if ~isempty(sensors(this, 'MEG'))
    timelock.grad = sensors(this, 'MEG');
end

if ~isempty(sensors(this, 'EEG'))
    timelock.elec = sensors(this, 'EEG');
end
