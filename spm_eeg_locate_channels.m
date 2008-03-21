function [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
% function [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
%
% Locates channels and generates mask for converting M/EEG data to nifti
% format ('analysis at sensor level')
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_locate_channels.m 1237 2008-03-21 14:54:07Z stefan $

% find channel positions on 2D plane
Cel = D.coor2D(D.meegchannels);

% check limits and stretch, if possible
dx = max(Cel(1,:)) - min(Cel(1,:));
dy = max(Cel(2,:)) - min(Cel(2,:));

if dx > 1 || dy > 1
    error('Coordinates not between 0 and 1');
end

scale = (1 - 10^(-6))/max(dx, dy);
Cel(1,:) = n*scale*(Cel(1,:) - min(Cel(1,:)) + eps) + 0.5;
Cel(2,:) = n*scale*(Cel(2,:) - min(Cel(2,:)) + eps) + 0.5;

% shift to middle
dx = n+0.5 -n*eps - max(Cel(1,:));
dy = n+0.5 -n*eps - max(Cel(2,:));
Cel(1,:) = Cel(1,:) + dx/2;
Cel(2,:) = Cel(2,:) + dy/2;

Cel = round(Cel)';

Cind = setdiff(D.meegchannels, D.badchannels);

[x, y] = meshgrid(1:n, 1:n);

if interpolate_bad
    % keep bad electrode positions in
    ch = convhull(Cel(:, 1), Cel(:, 2));
    Ic = find(inpolygon(x, y, Cel(ch, 1), Cel(ch, 2)));
else
    % or don't
    ch = convhull(Cel(Cind, 1), Cel(Cind, 2));
    Ic = find(inpolygon(x, y, Cel(Cind(ch), 1), Cel(Cind(ch), 2)));
end

Cel = Cel(Cind, :);

x = x(Ic); y = y(Ic);


