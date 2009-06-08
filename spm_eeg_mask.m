function spm_eeg_mask(S)
% Create a mask image for scalp-level contrasts based on an unsmoothed
% original image.
% FORMAT spm_eeg_mask(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%    image          - file name of an image containing an unsmoothed 
%                     M/EEG data in voxel-space
%    window         - start and end of a window in peri-stimulus time [ms]
%    outfile        - output file name
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_mask.m 3189 2009-06-08 17:09:31Z vladimir $

SVNrev = '$Rev: 3189 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG mask generation',0);

if nargin == 0
    S = [];
end


%-Input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'window')
    S.window = spm_input('start and end of window [ms or Hz]', '+1', 'r', '', 2);
end

if ~isfield(S, 'image')
    S.image = spm_select(1, 'image', 'Select M/EEG images (in voxel-space)');
end

if ~isfield(S, 'outfile')
    [filename, pathname] = uiputfile({'*.img;*.nii'}, 'Select the mask file');
    S.outfile = fullfile(pathname, filename);
end

spm('Pointer', 'Watch');


V = spm_vol(S.image);
Y = spm_read_vols(V);
Y = ~isnan(Y) & (Y~=0);

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

Y(: , :, 1:begsample)   = 0;
Y(: , :, endsample:end) = 0;

V.fname = S.outfile;

spm_write_vol(V, Y);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG mask generation: done'); spm('Pointer','Arrow');
