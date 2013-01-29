function spm_eeg_mask(S)
% Create a mask image for scalp-level contrasts.
% FORMAT spm_eeg_mask(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%    image          - file name of an image containing an unsmoothed 
%                     M/EEG data in voxel-space
%    timewin        - start and end of a window in peri-stimulus time [ms]
%    outfile        - output file name
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_mask.m 5219 2013-01-29 17:07:07Z spm $

SVNrev = '$Rev: 5219 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG mask generation'); spm('Pointer','Watch');

if nargin == 0
    S = [];
end

%-Input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'image')
    [S.image, sts] = spm_select(1, 'image', 'Select an unsmoothed M/EEG image (in voxel-space)');
    if ~sts, return; end
end

if ~isfield(S, 'timewin')
    S.timewin = spm_input('start and end of window [ms or Hz]', '+1', 'r', '', 2);
end

if ~isfield(S, 'outfile')
    [filename, pathname] = uiputfile({'*.img;*.nii'}, 'Select the mask file');
    [p, n, e] = fileparts(filename);
    if isempty(e), filename = [filename spm_file_ext]; end
    S.outfile = fullfile(pathname, filename);
end

V = spm_vol(S.image);
Y = spm_read_vols(V);
Y = ~isnan(Y) & (Y~=0);

Nt = size(Y, 3);

begsample =V.mat\[0 0 S.timewin(1) 1]';
begsample = begsample(3);

endsample = V.mat\[0 0 S.timewin(2) 1]';
endsample = endsample(3);

if any([begsample endsample] < 0) || ...
        any([begsample endsample] > Nt)
    error('The window is out of limits for the image.');
end

[junk,begsample] = min(abs(begsample-[1:Nt]));
[junk,endsample] = min(abs(endsample-[1:Nt]));

if begsample > 1
    Y(: , :, 1:(begsample-1))   = 0;
end

if endsample<size(Y, 3)
    Y(: , :, (endsample+1):end) = 0;
end

V.fname = S.outfile;

spm_write_vol(V, Y);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG mask generation: done'); spm('Pointer','Arrow');
