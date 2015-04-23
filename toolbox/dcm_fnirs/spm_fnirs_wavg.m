function wy = spm_fnirs_wavg(y, ons, dur)
% Average data across trials 
% FORMAT [wy] = spm_fnirs_wavg(y, fs, ons, dur)
%
% wy     time series averaged across trials
% 
% y       data (eg, optical density changes) 
% ons    onset of average window (eg, onset of tasks)
% dur    window size 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fnirs_wavg.m 6418 2015-04-23 08:40:08Z sungho $

n = length(ons); 
wy = zeros(dur, size(y,2)); 
for i = 1:n 
    wy = wy + y(ons(i):ons(i)+dur-1,:);
end
wy = wy./n; 
