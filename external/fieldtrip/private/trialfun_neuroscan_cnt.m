function [trl] = trialfun_ns_cnt(cfg);

% For continuous NeuroScan data, the trialdef structure should contain
%   cfg.trialdef.trigger  = number or list with triggers
%   cfg.trialdef.prestim  = pre-stimulus in seconds
%   cfg.trialdef.poststim = post-stimulus in seconds
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

% read the header and the event table from the continuous file
tmp = read_ns_cnt(cfg.datafile, 'ldheaderonly', 1);
fprintf('found %d events in continuous neuroscan file\n', length(tmp.event.frame));
eventindx = find(ismember(tmp.event.stimtype, cfg.trialdef.trigger));
trl(:,1) = tmp.event.frame(eventindx) - round(cfg.trialdef.prestim*tmp.rate) + 1;  % begin sample
trl(:,2) = tmp.event.frame(eventindx) + round(cfg.trialdef.poststim*tmp.rate) + 1; % end sample
trl(:,3) = -round(cfg.trialdef.prestim*tmp.rate);                                  % offset
fprintf('selected %d events based on triggercode\n', size(trl,1));
