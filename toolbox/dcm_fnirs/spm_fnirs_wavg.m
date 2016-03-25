function wy = spm_fnirs_wavg(y, ons, dur)
% Average data across trials 
% FORMAT [wy] = spm_fnirs_wavg(y, ons, dur)
%
% y       data (eg, optical density changes) 
% ons    onset of average window (eg, onset of tasks)
% dur    window size 
%
% wy     time series averaged across trials
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fnirs_wavg.m 6754 2016-03-25 06:44:58Z will $

n = length(ons); 
ns = size(y, 1); nch = size(y, 2); 
wy = NaN(dur, n, nch); 
for i = 1:n 
    eindx = ons(i)+dur -1; 
    if eindx < ns 
        wy(:, i, :) = y(ons(i):eindx, :); 
    else 
        ns_w = ns - ons(i) + 1; 
        wy(1:ns_w, i, :) = y(ons(i):ns, :); 
    end
end
wy = squeeze(nanmean(wy, 2)); 
