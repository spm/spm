function [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
% function [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
%
% Locates channels and generates mask for converting EEG data to analyze
% format on the scalp
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_locate_channels.m 776 2007-03-27 09:40:10Z stefan $

% load channel template file (contains location of channels)
Ctf = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

% find channel positions on 2D plane
Cel = Ctf.Cpos(:, D.channels.order(D.channels.eeg));

% check limits and stretch, if possible
dx = max(Cel(1,:)) - min(Cel(1,:));
dy = max(Cel(2,:)) - min(Cel(2,:));

if dx > 1 || dy > 1
    error('Check your channel template file: Coordinates not between 0 and 1');
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

Bad = [];
if isfield(D.channels, 'Bad')
    Bad = D.channels.Bad(:);
end

% For mapping indices
Itmp = zeros(1,length(D.channels.order));
Itmp(D.channels.eeg) = 1:length(D.channels.eeg);

Cind = setdiff(D.channels.eeg, Bad);

[x, y] = meshgrid(1:n, 1:n);
if interpolate_bad
    % keep bad electrode positions in
    ch = convhull(Cel(:, 1), Cel(:, 2));
    Ic = find(inpolygon(x, y, Cel(ch, 1), Cel(ch, 2)));
else
    % or don't
    ch = convhull(Cel(Itmp(Cind), 1), Cel(Itmp(Cind), 2));
    Ic = find(inpolygon(x, y, Cel(Itmp(Cind(ch)), 1), Cel(Itmp(Cind(ch)), 2)));
end

Cel = Cel(Itmp(Cind), :);

x = x(Ic); y = y(Ic);


