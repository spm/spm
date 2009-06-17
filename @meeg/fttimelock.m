function timelock = fttimelock(this)
% Method for converting meeg object to Fieldtrip timelock/freq struct
% FORMAT  timelock = fttimelock(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fttimelock.m 3210 2009-06-17 13:46:25Z vladimir $

timelock             = [];
timelock.label       = chanlabels(this);
timelock.label       = timelock.label(:);
if isequal(transformtype(this), 'time')
    timelock.dimord  = 'rpt_chan_time';
    timelock.fsample = fsample(this);
    timelock.trial   = permute(this.data.y(:, :, :), [3 1 2]);
    timelock.dimord  = 'rpt_chan_time';
elseif strncmpi(transformtype(this),'TF',2)
    timelock.dimord  = 'rpt_chan_freq_time';
    timelock.powspctrm = permute(this.data.y(:, :, :, :), [4 1 2 3]);
    timelock.freq    =  frequencies(this);
else
  error('Unknown transform type.');
end

timelock.time       = time(this);