function spm_eeg_firstlevel(S)
% function used for computing within-peristimulues time averages (contrasts)
% of M/EEG data in voxel-space
% FORMAT spm_eeg_firstlevel(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file defining the M/EEG parameters
% contrast1st   - struct with various entries:
%    images     - list of file names containing M/EEG data in voxel-space
%    window     - start and end of a window in peri-stimulus time [ms]
%    Pout       - output directory
%_______________________________________________________________________
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_firstlevel.m 2418 2008-10-30 10:33:39Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','M/EEG 1st level contrast setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

try
    contrast1st.window = S.contrast1st.window;
catch
    contrast1st.window = spm_input('start(s) and end(s) of window(s) [ms]', '+1', 'r', '', [inf 2]);
end


try
    contrast1st.fnames = S.contrast1st.fnames;
catch
    contrast1st.fnames = spm_select(inf, 'image', 'Select M/EEG images (in voxel-space)');
end

try
    contrast1st.Pout = S.contrast1st.Pout;
catch
    contrast1st.Pout = uigetdir('.', 'Select output directory');
end

% check input
w = contrast1st.window;
Nc = size(w, 1);

if any(w(:, 1) < time(D, 1, 'ms'))
    error('Start of time window must be later than %d ms.', time(D, 1));
    return;
end

if any(w(:, 2) > time(D, D.nsamples, 'ms'))
    error('End of time window must be earlier than %d ms.', time(D, D.nsamples));
    return;
end

if any(w(:, 1) > w(:, 2))
    error('Start of time window must be earlier than its end.');
    return;
end

spm('Pointer', 'Watch'); drawnow;

% compute contrasts
spm('defaults','EEG');
global defaults
defaults.analyze.flip = 0;

C = zeros(D.nsamples, Nc);

for i = 1:Nc
    tsample(i, 1) = indsample(D, w(i, 1)/1000);
    tsample(i, 2) = indsample(D, w(i, 2)/1000);
    C(tsample(1):tsample(2), i) = 1./(tsample(i, 2) - tsample(i, 1) + 1);
end

fnames = cellstr(contrast1st.fnames);

for j = 1:length(fnames) % over files

    % map file
    Vbeta = nifti(deblank(fnames{j}));

    % cd to target directory
    cd(contrast1st.Pout);

    for i = 1:Nc % over contrasts

        % code taken from spm_contrasts
        fprintf('\t%-32s: %-10s%20s', sprintf('file %s, contrast %d', fnames{j}, i),...
            '(spm_add)','...initialising') %-#

        %-Prepare handle for contrast image
        %-----------------------------------------------------------
        descrip = sprintf('SPM contrast - average from %d to %d ms',...
                w(i, 1), w(i, 2));

        % prepare nifti image (the usual spm_add doesn't seem to work for
        % many input files under windows)
        Vcon = nifti;
        Vcon.descrip = descrip;

        Dcon = file_array;
        Dcon.fname = sprintf('%s_con_%04d.img', spm_str_manip(fnames{j}, 'rt'), i);
        Dcon.dtype = spm_type('float32');
        Dcon.offset  = ceil(348/8)*8;
        Dcon.dim = Vbeta.dat.dim(1:2);

        %-Write image
        %-----------------------------------------------------------
        fprintf('%s%20s', repmat(sprintf('\b'),1,20),'...computing')%-#

        d = zeros(Vbeta.dat.dim(1:2));
        for k = 1:Vbeta.dat.dim(3)
            d = d + Vbeta.dat(:, : ,k)*C(k,i);
        end

        Dcon(:,:) = d;

        Vcon.dat = Dcon;
        create(Vcon);

        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(Vcon.dat.fname,'t')))%-#


    end
end

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');
