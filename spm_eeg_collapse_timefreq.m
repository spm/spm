function images = spm_eeg_collapse_timefreq(S)
% Compute within-peristimulus time (or frequency) averages (contrasts) of M/EEG data in voxel-space
% FORMAT images = spm_eeg_collapse_timefreq(S)
%
% S         - input structure 
% fields of S:
%    images  - list of file names containing M/EEG data in voxel-space
%    timewin - C x 2 matrix of start(s) and end(s) of a window in peri-stimulus 
%              time [ms] (or frequency (Hz))
%    prefix  - prefix for the averaged images
%_____________________________________________________________________________________________
% Copyright (C) 2006-2013 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_collapse_timefreq.m 5194 2013-01-18 15:04:19Z vladimir $

SVNrev = '$Rev: 5194 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG Collapse time',0);

if ~isfield(S, 'prefix'),       S.prefix   = 'l';           end
if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end

spm('Pointer', 'Watch');

%-Compute contrasts
%--------------------------------------------------------------------------
Nf = size(S.images, 1);
Nc = size(S.timewin, 1);

if ischar(S.images)
    fnames = cellstr(S.images);
else
    fnames = S.images;
end

ind1   = find(S.timewin(:, 1) == -Inf);
S.timewin(ind1, 1)   = 0;

indend = find(S.timewin(:, 2) ==  Inf);
S.timewin(indend, 2) = 0;

images = {};

spm_progress_bar('Init', Nf, 'Averaging in images');
if Nf > 100, Ibar = floor(linspace(1, Nf, 100));
else Ibar = 1:Nf; end

for j = 1:Nf % over files

    Vbeta = nifti(fnames{j});

    Nt = size(Vbeta.dat, 3); % Number of time frames

    begsample = Vbeta.mat\[zeros(2, Nc); S.timewin(:, 1)'; ones(1, Nc)];
    begsample = begsample(3, :);
    begsample(ind1) = 1;

    endsample = Vbeta.mat\[zeros(2, Nc); S.timewin(:, 2)'; ones(1, Nc)];
    endsample = endsample(3, :);
    endsample(indend) = Nt;

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
            spm_file(fnames{j}, 'basename'), i), '...initialising');    %-#

        %-Write contrast image header
        %------------------------------------------------------------------
        Vcon               = Vbeta;
        Vcon.mat(3,3:4)    = [1.0 0.0];
        Vcon.mat0          = Vcon.mat;
        Vcon.dat.fname     = spm_file(fnames{j}, 'basename', sprintf('%s%s_con_%04d', S.prefix, spm_file(fnames{j},'basename'), i), 'ext', spm_file_ext);
        Vcon.dat.scl_slope = 1.0;
        Vcon.dat.scl_inter = 0.0;
        Vcon.dat.dtype     = 'float32-le';
        Vcon.dat.offset    = 0;
        Vcon.dat.dim       = Vbeta.dat.dim(1:2);
        Vcon.descrip       = sprintf('SPM contrast - average from %d to %d ms',...
            S.timewin(i, 1), S.timewin(i, 2));
        create(Vcon);
        
        images{end+1}      = Vcon.dat.fname;

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

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG Collapse time: done'); spm('Pointer','Arrow');
