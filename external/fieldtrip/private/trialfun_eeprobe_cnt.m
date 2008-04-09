function [trl] = trialfun_eeprobe_cnt(cfg);

% For continuous EEProbe data, the trialdef structure should contain
%   cfg.trialdef.trigger  = 2             trigger code, can be a list
%   cfg.trialdef.prestim  = 0.300         latency in seconds
%   cfg.trialdef.poststim = 0.700         latency in seconds
%   cfg.trialdef.trgfile  = filename of the trigger file
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

% read external trigger file in EEP format
trg = read_eep_trg(cfg.trialdef.trgfile);
trl = [];
for i=1:length(trg)
  if ismember(trg(i).type,cfg.trialdef.trigger)
    trlbegin  = round((trg(i).time/1000-cfg.trialdef.prestim)*hdr.Fs)+1;
    trlend    = round((trg(i).time/1000+cfg.trialdef.poststim)*hdr.Fs)+1;
    trloffset = round(-cfg.trialdef.prestim*hdr.Fs);
    trl = [trl; trlbegin trlend trloffset];
  end
end
fprintf('%d triggers converted into %d trials\n', length(trg), size(trl,1));
