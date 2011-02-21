function spm_eeg_img2maps(S)
% Make a series of scalp maps from data in an image
% FORMAT  spm_eeg_img2maps(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%    D              - M/EEG dataset containing the sensor locations
%    image          - file name of an image containing M/EEG data in voxel-space
%    window         - start and end of a window in peri-stimulus time [ms]
%    clim           - color limits of the plot
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_img2maps.m 4204 2011-02-21 12:36:08Z vladimir $

SVNrev = '$Rev: 4204 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Plot scalp maps');

if nargin == 0
    S = [];
end

%% -Input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'image')
    [S.image, sts] = spm_select(1, 'image', 'Select an M/EEG image (in voxel-space)');
    if ~sts, return; end
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~isfield(S, 'window')
    S.window = spm_input('start and end of window [ms or Hz]', '+1', 'r', '', 2);
end

V = spm_vol(S.image);
Y = spm_read_vols(V);

Nt = size(Y, 3);

begsample = inv(V.mat)*[0 0 S.window(1) 1]';
begsample = begsample(3);

endsample = inv(V.mat)*[0 0 S.window(2) 1]';
endsample = endsample(3);

if any([begsample endsample] < 0) || ...
        any([begsample endsample] > Nt)
    error('The window is out of limits for the image.');
end

[junk begsample] = min(abs(begsample-[1:Nt]));
[junk endsample] = min(abs(endsample-[1:Nt]));

Y = squeeze(mean(Y(: , :, begsample:endsample), 3));

if any(diff(size(Y)))
    error('The image should be square');
else
    n = size(Y, 1);
end

%-Get channel indices and coordinates
%--------------------------------------------------------------------------
[Cel, Cind] = spm_eeg_locate_channels(D, n, 1);

modality    = spm_eeg_modality_ui(D, 1, 1);

goodchan    = find(ismember(Cind, setdiff(D.meegchannels(modality), D.badchannels)));

Cel         = Cel(goodchan, :);
Cind        = Cind(goodchan);

if isempty(Cind)
    error('No good channels to plot');
end

Y = Y(sub2ind(size(Y), Cel(:, 1), Cel(:, 2)));


%%
% SPM graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph)

if ~isfield(S, 'clim')
    in.max = max(abs(Y));
    in.min = -in.max;
else
    in.min = S.clim(1);
    in.max = S.clim(2);
end

in.type = modality;
in.f = Fgraph;
in.ParentAxes = axes;

spm_eeg_plotScalpData(Y, D.coor2D(Cind), D.chanlabels(Cind), in);
%%
%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Plot scalp maps: done');
