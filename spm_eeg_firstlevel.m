function spm_eeg_firstlevel(S)
% Compute within-peristimulus time averages (contrasts) of M/EEG data in voxel-space
% FORMAT spm_eeg_firstlevel(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%    images         - list of file names containing M/EEG data in voxel-space
%    window         - start and end of a window in peri-stimulus time [ms]
%    Pout           - output directory
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_firstlevel.m 3138 2009-05-20 14:32:53Z vladimir $

SVNrev = '$Rev: 3138 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG 1st level contrast setup',0);

if nargin == 0
    S = [];
end

%-Backward compatibility
%--------------------------------------------------------------------------
if isfield(S, 'contrast1st')
    S = S.contrast1st;
end

if isfield(S, 'fnames')
    S.images = S.fnames;
end

%-Input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'window')
    S.window = spm_input('start(s) and end(s) of window(s) [ms]', '+1', 'r', '', [Inf 2]);
end

if ~isfield(S, 'images')
    S.images = spm_select(Inf, 'image', 'Select M/EEG images (in voxel-space)');
end

if ~isfield(S, 'Pout')
    S.Pout = spm_select(1, 'dir', 'Select output directory');
end

spm('Pointer', 'Watch');

%-Change to target directory
%--------------------------------------------------------------------------
swd = pwd;
cd(S.Pout);

%-Compute contrasts
%--------------------------------------------------------------------------
Nf = size(S.images, 1);
Nc = size(S.window, 1);

fnames = cellstr(S.images);

spm_progress_bar('Init', Nf, 'First level contrasts');
if Nf > 100, Ibar = floor(linspace(1, Nf, 100));
else Ibar = 1:Nf; end

for j = 1:Nf % over files

    Vbeta = nifti(fnames{j});

    Nt = size(Vbeta.dat, 3); % Number of time frames

    begsample = inv(Vbeta.mat)*[zeros(2, Nc); S.window(:, 1)'; ones(1, Nc)];
    begsample = begsample(3, :);

    endsample = inv(Vbeta.mat)*[zeros(2, Nc); S.window(:, 2)'; ones(1, Nc)];
    endsample = endsample(3, :);

    if any([begsample endsample] < 0) || ...
            any([begsample endsample] > Nt)
        error(['The window is out of limits for image ' fnames{j}]);
    end

    for i = 1:Nc % over contrasts
        C  = zeros(Nt, 1);

        tsample = [];
        [junk, tsample(1)] = min(abs((1:Nt) - begsample(i)));
        [junk, tsample(2)] = min(abs((1:Nt) - endsample(i)));

        C(tsample(1):tsample(2)) = 1./(tsample(2) - tsample(1) + 1);


        fprintf('%-40s: %30s', sprintf('file %s, contrast %d', ...
            spm_str_manip(fnames{j}, 'rt'), i), '...initialising');     %-#

        %-Write contrast image header
        %------------------------------------------------------------------
        Vcon               = Vbeta;
        Vcon.mat(3,3:4)    = [1.0 0.0];
        Vcon.mat0          = Vcon.mat;
        Vcon.dat.fname     = sprintf('%s_con_%04d.img', spm_str_manip(fnames{j}, 'rt'), i);
        Vcon.dat.scl_slope = 1.0;
        Vcon.dat.scl_inter = 0.0;
        Vcon.dat.dtype     = 'float32-le';
        Vcon.dat.offset    = 0;
        Vcon.dat.dim       = Vbeta.dat.dim(1:2);
        Vcon.descrip       = sprintf('SPM contrast - average from %d to %d ms',...
            S.window(i, 1), S.window(i, 2));
        create(Vcon);

        %-Compute contrast
        %------------------------------------------------------------------
        fprintf('%s%30s', repmat(sprintf('\b'),1,30),'...computing');   %-#

        d = zeros(Vbeta.dat.dim(1:2));
        for k = 1:Vbeta.dat.dim(3)
            d = d + Vbeta.dat(:, : ,k) * C(k);
        end

        %-Write contrast image
        %------------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...writing');      %-#

        Vcon.dat(:,:) = d;

        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...written');    %-#

    end

    if ismember(j, Ibar), spm_progress_bar('Set', j); end

end

cd(swd);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG 1st level contrast: done'); spm('Pointer','Arrow');
