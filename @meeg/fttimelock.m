function timelock = fttimelock(this, chanind, timeind, trialind, freqind)
% Method for converting meeg object to Fieldtrip timelock/freq struct
% FORMAT timelock = fttimelock(this, chanind, timeind, trialind, freqind)
%
% The method support both time and TF data and outputs different variants
% of timelock or freq FT struct depending on the dataset type and requested
% data dimensions.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ~islinked(this)
    error('There is no data linked to the object');
end

if nargin < 2 || ~isnumeric(chanind)
    chanind = 1:nchannels(this);
end

if nargin < 3 || ~isnumeric(timeind)
    timeind = 1:nsamples(this);
end

if nargin < 4 || ~isnumeric(trialind)
    trialind = 1:ntrials(this);
end

if strncmpi(transformtype(this),'TF',2) && ...
        (nargin < 5 || isempty(freqind))
     freqind = 1:nfrequencies(this);
end

timelock             = [];
timelock.label       = chanlabels(this, chanind)';

if isequal(transformtype(this), 'time')
   if isequal(type(this), 'continuous')
            error('For continuous data use ftraw method');
   end
   
   if isequal(type(this), 'single') || length(trialind)>1
       timelock.dimord  = 'rpt_chan_time';
       timelock.trial   =  permute(subsref(this, substruct('()', {chanind, timeind, trialind})), [3 1 2]);
   else
       timelock.dimord  = 'chan_time';
       timelock.avg     =  spm_squeeze(subsref(this, substruct('()', {chanind, timeind, trialind})), 3);
   end
   
   timelock.time       = time(this, timeind);
   
elseif strncmpi(transformtype(this),'TF',2)
    if length(timeind)>1
        if isequal(type(this), 'single') || length(trialind)>1
            timelock.dimord    = 'rpt_chan_freq_time';
            timelock.powspctrm = permute(subsref(this, substruct('()', {chanind, freqind, timeind, trialind})), [4 1 2 3]);
        else
            timelock.dimord     = 'chan_freq_time';
            timelock.powspctrm  =  spm_squeeze(subsref(this, substruct('()', {chanind, freqind, timeind, trialind})), 3);
        end
        
        timelock.time       = time(this, timeind);
    else
        if isequal(type(this), 'single') || length(trialind)>1
            timelock.dimord    = 'rpt_chan_freq';
            timelock.powspctrm = spm_squeeze(permute(subsref(this, substruct('()', {chanind, freqind, timeind, trialind})), [4 1 2 3]), 4);
        else
            timelock.dimord     = 'chan_freq';
            timelock.powspctrm  =  spm_squeeze(subsref(this, substruct('()', {chanind, freqind, timeind, trialind})), [3 4]);
        end
    end
    
    timelock.freq      = frequencies(this, freqind);
else
    error('Unknown transform type.');
end   

if length(trialind)>1
    
    clist      =  condlist(this);
    condlabels = conditions(this, trialind);
    timelock.trialinfo = 0*trialind;
    
    for k = 1:numel(clist)
        fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
        sel = strcmp(clist{k}, condlabels);
        timelock.trialinfo(sel) = k;
    end
    
end

if ~isempty(sensors(this, 'MEG'))
    timelock.grad = sensors(this, 'MEG');
end

if ~isempty(sensors(this, 'EEG'))
    timelock.elec = sensors(this, 'EEG');
end

hdr = [];
hdr.Fs          = fsample(this);
hdr.nChans      = length(chanind);
hdr.nSamples    = length(timeind);
hdr.nSamplesPre = sum(time(this, timeind)<0);
hdr.nTrials     = length(trialind);
hdr.label       = timelock.label;
hdr.chanunit    = units(this, chanind);

spmtype      = chantype(this, chanind);
dictionary   = spm_eeg_spmft_chan_dictionary;
[sel1, sel2] = spm_match_str(spmtype, dictionary(:, 2));
        
hdr.chantype = dictionary(sel2, 1)';

if isfield(this.other, 'origheader')
    hdr.orig = this.other.origheader;
end

timelock.hdr = hdr;
