function spm_eeg_firstlevel(S)
% Compute within-peristimulus time averages (contrasts) of M/EEG data in voxel-space
% FORMAT spm_eeg_firstlevel(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%   S.D      - MEEG object or filename of M/EEG mat-file
%   S.contrast1st   - struct with various entries:
%    images         - list of file names containing M/EEG data in voxel-space
%    window         - start and end of a window in peri-stimulus time [ms]
%    Pout           - output directory
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_firstlevel.m 2829 2009-03-05 12:05:07Z guillaume $

SVNrev = '$Rev: 2829 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG 1st level contrast setup',0);

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, return; end
end

D = spm_eeg_load(D);

%-Input parameters
%--------------------------------------------------------------------------
try
    contrast1st.window = S.contrast1st.window;
catch
    contrast1st.window = spm_input('start(s) and end(s) of window(s) [ms]', '+1', 'r', '', [Inf 2]);
end

try
    contrast1st.fnames = S.contrast1st.fnames;
catch
    contrast1st.fnames = spm_select(Inf, 'image', 'Select M/EEG images (in voxel-space)');
end

try
    contrast1st.Pout = S.contrast1st.Pout;
catch
    contrast1st.Pout = uigetdir(pwd, 'Select output directory');
end

%-Check input
%--------------------------------------------------------------------------
w  = contrast1st.window;

if any(w(:, 1) < time(D, 1, 'ms'))
    error('Start of time window must be later than %d ms.', time(D, 1));
end

if any(w(:, 2) > time(D, D.nsamples, 'ms'))
    error('End of time window must be earlier than %d ms.', time(D, D.nsamples));
end

if any(w(:, 1) > w(:, 2))
    error('Start of time window must be earlier than its end.');
end

spm('Pointer', 'Watch');

%-Change to target directory
%--------------------------------------------------------------------------
swd = pwd;
cd(contrast1st.Pout);

%-Compute contrasts
%--------------------------------------------------------------------------
Nc = size(w, 1);
C  = zeros(D.nsamples, Nc);

for i = 1:Nc
    tsample(1) = indsample(D, w(i, 1)/1000);
    tsample(2) = indsample(D, w(i, 2)/1000);
    C(tsample(1):tsample(2), i) = 1./(tsample(2) - tsample(1) + 1);
end

fnames = cellstr(contrast1st.fnames);

spm_progress_bar('Init', length(fnames), 'First level contrasts');
if length(fnames) > 100, Ibar = floor(linspace(1, length(fnames), 100));
else Ibar = 1:length(fnames); end

for j = 1:length(fnames) % over files

    Vbeta = nifti(fnames{j});

    for i = 1:Nc % over contrasts

        fprintf('%-40s: %30s', sprintf('file %s, contrast %d', ...
            spm_str_manip(fnames{j}, 'rt'), i), '...initialising');     %-#

        %-Write contrast image header
        %------------------------------------------------------------------
        Vcon               = Vbeta;
        Vcon.dat.fname     = sprintf('%s_con_%04d.img', spm_str_manip(fnames{j}, 'rt'), i);
        Vcon.dat.scl_slope = 1.0;
        Vcon.dat.scl_inter = 0.0;
        Vcon.dat.dtype     = 'float32-le';
        Vcon.dat.offset    = 0;
        Vcon.dat.dim       = Vbeta.dat.dim(1:2);
        Vcon.descrip       = sprintf('SPM contrast - average from %d to %d ms',...
                                w(i, 1), w(i, 2));
        create(Vcon);

        %-Compute contrast
        %------------------------------------------------------------------
        fprintf('%s%30s', repmat(sprintf('\b'),1,30),'...computing');   %-#

        d = zeros(Vbeta.dat.dim(1:2));
        for k = 1:Vbeta.dat.dim(3)
            d = d + Vbeta.dat(:, : ,k) * C(k,i);
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
