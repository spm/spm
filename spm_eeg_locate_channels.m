function [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
% Locate channels and generate mask for converting M/EEG data into images
% FORMAT [Cel, Cind, x, y] = spm_eeg_locate_channels(D, n, interpolate_bad)
%
% D               - M/EEG object
% n               - number of voxels in each direction
% interpolate_bad - flag (1/0), whether bad channels should be interpolated
%                   or masked out
%
% Cel             - coordinates of good channels in new coordinate system
% Cind            - the indices of these channels in the total channel
%                   vector
% x, y            - x and y coordinates which support data
%
%__________________________________________________________________________
%
% Locates channels and generates mask for converting M/EEG data to NIfTI
% format ('analysis at sensor level'). If flag interpolate_bad is set to 1,
% the returned x,y-coordinates will include bad sensor position. If
% interpolate_bad is 0, these locations are masked out if the sensor are 
% located at the edge of the setup (where the data cannot be well 
% interpolated).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_locate_channels.m 2830 2009-03-05 17:27:34Z guillaume $

% put into n x n grid
%--------------------------------------------------------------------------
[x, y] = meshgrid(1:n, 1:n);

if interpolate_bad
    % keep bad electrode positions in
    %----------------------------------------------------------------------
    Cel = scale_coor(D.coor2D(D.meegchannels), n);
else
    % bad electrodes are masked out if located at the edge of the setup
    %----------------------------------------------------------------------
    Cel = scale_coor(D.coor2D(setdiff(D.meegchannels, D.badchannels)), n);
end

ch = convhull(Cel(:, 1), Cel(:, 2));
Ic = find(inpolygon(x, y, Cel(ch, 1), Cel(ch, 2)));

Cel = scale_coor(D.coor2D(setdiff(D.meegchannels, D.badchannels)), n);
Cind = setdiff(D.meegchannels, D.badchannels);

x = x(Ic); y = y(Ic);

%==========================================================================
% scale_coor
%==========================================================================
function Cel = scale_coor(Cel, n)

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

% 2D coordinates in voxel-space (incl. badchannels)
Cel = round(Cel)';
