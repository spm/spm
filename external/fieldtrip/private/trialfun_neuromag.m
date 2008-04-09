function [trl] = trialfun_neuromag(cfg);

% For Neuromag data, the trialdef structure should contain
%   cfg.trialdef.trgchan  = channel label, e.g. 'STI 001'
%   cfg.trialdef.prestim  = 0.300         latency in seconds
%   cfg.trialdef.poststim = 0.700         latency in seconds
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

chanindx = match_str(hdr.label, cfg.trialdef.trgchan);
nsamples = hdr.nTrials * hdr.nSamples;
if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end
read_fcdc_data(cfg.datafile, hdr, 1, nsamples, chanindx, iscontinuous);
trigindx = find(trigger & [0 diff(trigger)]);
trl(:,1) = trigindx(:) - round(cfg.trialdef.prestim*hdr.Fs);  % define begin of each trial (in samples)
trl(:,2) = trigindx(:) + round(cfg.trialdef.poststim*hdr.Fs); % define end of each trial (in samples)
trl(:,3) = round(-cfg.trialdef.prestim*hdr.Fs);               % define the trial offset relative to latency zero (in samples)
trl(find(trl(:,1)<1), :) = [];                                % remove trials beginning before the start of the file
trl(find(trl(:,2)>hdr.nTrials*hdr.nSamples), :) = [];         % remove trials after the end of the file
fprintf('%d triggers converted into %d trials\n', length(trigindx), size(trl,1));

