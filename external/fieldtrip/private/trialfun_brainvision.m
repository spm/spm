function [trl] = trialfun_brainvision(cfg);

% For continuous BrainVision data, the trialdef structure should contain
%   cfg.trialdef.stim     = 2             stimulus trigger code, can be a list
%   cfg.trialdef.resp     = 2             response trigger code, can be a list
%   cfg.trialdef.segment  = 'no' / 'yes'  use the segment markers as onset
%   cfg.trialdef.timezero = 'no' / 'yes'  use the "Time 0" markers as onset
%   cfg.trialdef.prestim  = 0.300         pre-marker latency in seconds
%   cfg.trialdef.poststim = 0.700         post-marker latency in seconds
%   cfg.trialdef.trgfile  = filename of the marker file
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

% read the external marker file in BrainVision format
[stim, resp, segment, timezero] = read_brainvision_vmrk(cfg.trialdef.trgfile);
% set the default marker elements that can be converted into trials
if ~isfield(cfg.trialdef, 'stim')
  cfg.trialdef.stim = [];
end
if ~isfield(cfg.trialdef, 'resp')
  cfg.trialdef.resp = [];
end
if ~isfield(cfg.trialdef, 'segment')
  cfg.trialdef.segment = 'no';
end
if ~isfield(cfg.trialdef, 'timezero')
  cfg.trialdef.timezero = 'no';
end
trl = [];
for i=1:size(stim,1)
  if ismember(stim(i,1),cfg.trialdef.stim)
    trlbegin  = round((stim(i,2)/1000-cfg.trialdef.prestim)*hdr.Fs)+1;
    trlend    = round((stim(i,2)/1000+cfg.trialdef.poststim)*hdr.Fs)+1;
    trloffset = round(-cfg.trialdef.prestim*hdr.Fs);
    trl = [trl; trlbegin trlend trloffset];
  end
end
fprintf('%d stimulus markers converted into %d trials\n', size(stim,1), size(trl,1));
tmp = size(trl,1);
for i=1:size(resp,1)
  if ismember(resp(i,1),cfg.trialdef.resp)
    trlbegin  = round((resp(i,2)/1000-cfg.trialdef.prestim)*hdr.Fs)+1;
    trlend    = round((resp(i,2)/1000+cfg.trialdef.poststim)*hdr.Fs)+1;
    trloffset = round(-cfg.trialdef.prestim*hdr.Fs);
    trl = [trl; trlbegin trlend trloffset];
  end
end
fprintf('%d response markers converted into %d trials\n', size(resp,1), size(trl,1)-tmp);
tmp = size(trl,1);
if strcmp(cfg.trialdef.segment, 'yes')
  for i=1:size(segment,1)
    trlbegin  = round((segment(i)/1000-cfg.trialdef.prestim)*hdr.Fs)+1;
    trlend    = round((segment(i)/1000+cfg.trialdef.poststim)*hdr.Fs)+1;
    trloffset = round(-cfg.trialdef.prestim*hdr.Fs);
    trl = [trl; trlbegin trlend trloffset];
  end
end
fprintf('%d segment markers converted into %d trials\n', size(segment,1), size(trl,1)-tmp);
tmp = size(trl,1);
if strcmp(cfg.trialdef.timezero, 'yes')
  for i=1:size(timezero,1)
    trlbegin  = round((timezero(i)/1000-cfg.trialdef.prestim)*hdr.Fs)+1;
    trlend    = round((timezero(i)/1000+cfg.trialdef.poststim)*hdr.Fs)+1;
    trloffset = round(-cfg.trialdef.prestim*hdr.Fs);
    trl = [trl; trlbegin trlend trloffset];
  end
end
fprintf('%d "Time 0" markers converted into %d trials\n', size(timezero,1), size(trl,1)-tmp);
