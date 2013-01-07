function [Cel, x, y] = spm_eeg_locate_channels(D, n, Cind)
% Locate channels and generate mask for converting M/EEG data into images
% FORMAT [Cel, x, y] = spm_eeg_locate_channels(D, n, channels)
%
% D               - M/EEG object
% n               - number of voxels in each direction
% Cind            - the indices of channels in the total channel
%                   vector
%
% Cel             - coordinates of channels in new coordinate system
% x, y            - x and y coordinates which support data
%
%__________________________________________________________________________
%
% Locates channels and generates mask for converting M/EEG data to NIfTI
% format ('analysis at sensor level'). 
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_locate_channels.m 5177 2013-01-07 11:36:08Z vladimir $

% put into n x n grid
%--------------------------------------------------------------------------
[x, y] = meshgrid(1:n, 1:n);

Cel  = scale_coor(coor2D(D, Cind), n);

ch = convhull(Cel(:, 1), Cel(:, 2));
Ic = find(inpolygon(x, y, Cel(ch, 1), Cel(ch, 2)));

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
