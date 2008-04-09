function [trl] = trialfun_ns_eeg(cfg);

% For epoched NeuroScan data, the trialdef structure should not contain anything
%   cfg.trialdef = []
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

for k=1:hdr.nsweeps
  trl(k,1) = 1+(k-1)*hdr.npnt;
  trl(k,2) = k*hdr.npnt;
  trl(k,3) = (hdr.rate*hdr.xmin/1000);
end
fprintf('found %d trials in epoched neuroscan file\n', size(trl,1));
fprintf('looking for manually accepted/rejected trials ');
accept = ones(1,hdr.nsweeps);
for k=1:hdr.nsweeps
  tmp = read_ns_eeg(cfg.datafile, k);
  accept(k) = tmp.sweep.accept;
  if accept(k)
    fprintf('.');
  else
    fprintf('R');
  end
end
fprintf('\n');
trl = trl(find(accept), :);
fprintf('accepted %d trials, rejected %d trials\n', sum(accept==1), sum(accept==0));
