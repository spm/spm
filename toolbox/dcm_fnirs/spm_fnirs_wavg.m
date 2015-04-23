function wy = spm_fnirs_wavg(y, ons, dur)
% Average data across trials
% FORMAT [wy] = spm_fnirs_wavg(y, ons, dur)
%
% y    - data (eg, optical density changes)
% ons  - onset of average window (eg, onset of tasks)
% dur  - window size
%
% wy   - time series averaged across trials
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fnirs_wavg.m 6422 2015-04-23 16:55:52Z spm $

n = length(ons);
wy = zeros(dur, size(y,2));
for i = 1:n
    wy = wy + y(ons(i):ons(i)+dur-1,:);
end
wy = wy./n;
