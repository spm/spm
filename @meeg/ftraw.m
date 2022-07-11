function raw = ftraw(this, chanind, timeind, trialind)
% Method for converting meeg object to Fieldtrip raw struct
% FORMAT raw  =  ftraw(this, chanind, timeind, trialind)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ~islinked(this)
    error('There is no data linked to the object');
end

if ~isequal(transformtype(this), 'time')
    raw = fttimelock(this, chanind, timeind, trialind);
    return;
end

% chanind == 0 is accepted for backward compatibility
if nargin < 2 || ~isnumeric(chanind) || isequal(chanind, 0)
    chanind = 1:nchannels(this);
end

if nargin < 3 || ~isnumeric(timeind)
    timeind = 1:nsamples(this);
end

if nargin < 4 || ~isnumeric(trialind)
    trialind = 1:ntrials(this);
end

raw = [];

raw.label   = chanlabels(this, chanind)';

raw.trial   = cell(1, length(trialind));

for i =  1:length(trialind) 
    raw.trial{i} = subsref(this, substruct('()', {chanind, timeind, trialind(i)}));
end

raw.time = repmat({time(this, timeind)}, 1, length(trialind));

clist      =  condlist(this);

condlabels = conditions(this, trialind);

raw.trialinfo = 0*trialind;
for k = 1:numel(clist)
  fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
  sel = strcmp(clist{k}, condlabels);
  raw.trialinfo(sel) = k;
end

if ~isempty(sensors(this, 'MEG'))
    raw.grad = sensors(this, 'MEG');
end

if ~isempty(sensors(this, 'EEG'))
    raw.elec = sensors(this, 'EEG');
end

onsets = trialonset(this, trialind);

if all(onsets>0)
    onsets = round(onsets(:)*fsample(this));
    raw.sampleinfo = [onsets+timeind(1) onsets+timeind(end)]-1;
end

hdr = [];
hdr.Fs          = fsample(this);
hdr.nChans      = length(chanind);
hdr.nSamples    = length(timeind);
hdr.nSamplesPre = sum(time(this, timeind)<0);
hdr.nTrials     = length(trialind);
hdr.label       = raw.label;
hdr.chanunit    = units(this, chanind);

spmtype      = chantype(this, chanind);
dictionary   = spm_eeg_spmft_chan_dictionary;
[sel1, sel2] = spm_match_str(spmtype, dictionary(:, 2));
        
hdr.chantype = dictionary(sel2, 1)';

if isfield(this.other, 'origheader')
    hdr.orig = this.other.origheader;
end

raw.hdr = hdr;

