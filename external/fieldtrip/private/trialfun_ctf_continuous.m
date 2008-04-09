function [trl] = trialfun_ctf_continuous(cfg);

% For continuous CTF data, the trialdef structure should contain
%   cfg.trialdef.trigger  = 2             trigger code, can be an array
%   cfg.trialdef.prestim  = 0.300         latency in seconds
%   cfg.trialdef.poststim = 0.700         latency in seconds
%   cfg.trialdef.panel    = 'front' or 'back' (optional, back is the default)

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

% read the triggers from the STIM channel
fprintf('reading stimulus channel and detecting trigger onset\n');
[backpanel, frontpanel] = read_ctf_trigger(cfg.dataset);
if isfield(cfg.trialdef, 'panel')
  if strcmp(cfg.trialdef.panel, 'back')
    trigger = backpanel;
  elseif strcmp(cfg.trialdef.panel, 'front')
    trigger = frontpanel;
  elseif strcmp(cfg.trialdef.panel, 'both')
    % recombine the triggers on back and frontpanel into a single trigger value
    trigger = frontpanel + (2^16)*backpanel;
  else
    error('unrecognized option for cfg.trialdef.panel');
  end
else
  % default is to use the backpanel where the Presentations PC is plugged into
  trigger = backpanel;
end
fprintf('total number of triggers is %d\n', length(find(trigger)));
% find the samples at which the specified trigger occured
trigindx = find(ismember(trigger, cfg.trialdef.trigger));
% convert the trigger moments into trials
trl(:,1) = trigindx(:) - round(cfg.trialdef.prestim*hdr.Fs);  % define begin of each trial (in samples)
trl(:,2) = trigindx(:) + round(cfg.trialdef.poststim*hdr.Fs); % define end of each trial (in samples)
trl(:,3) = round(-cfg.trialdef.prestim*hdr.Fs);               % define the trial offset relative to latency zero (in samples)
trl(find(trl(:,1)<1), :) = [];                                % remove trials beginning before the start of the file
trl(find(trl(:,2)>hdr.nTrials*hdr.nSamples), :) = [];         % remove trials after the end of the file
fprintf('%d triggers converted into %d trials\n', length(trigindx), size(trl,1));
