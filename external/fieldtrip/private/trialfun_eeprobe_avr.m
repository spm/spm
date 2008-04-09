function [trl] = trialfun_eeprobe_avr(cfg);

% For averaged EEProbe data, the trialdef structure should not contain anything
%   cfg.trialdef = []
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

trl(1) = 1;
trl(2) = hdr.nSamples;
trl(3) = -hdr.nSamplesPre;
