function D = spm_eeg_downsample(S)
% function used for down-sampling EEG/MEG data
% FORMAT D = spm_eeg_downsample(S)
% 
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file
% Radc_new  - new sampling rate
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 182 2005-05-27 17:44:15Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG downsample setup',0);

try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

try
    Radc_new = S.Radc_new;
catch
    str = 'New sampling rate';
    YPos = -1;
    while 1
        if YPos == -1
            YPos = '+1';
        end
        [Radc_new, YPos] = spm_input(str, YPos, 'r');
        if Radc_new < D.Radc, break, end
        str = sprintf('Sampling rate must be less than original (%d)', round(D.Radc));
    end
end

spm('Pointer', 'Watch');drawnow;

% Prepare for writing data
D.fnamedat = ['d' D.fnamedat];
fpd = fopen(fullfile(P, D.fnamedat), 'w');

% treat continuous and epoched data differently because of different
% scaling method

if size(D.data, 3) > 1
    % epoched
    D.scale.dim = [1 3];
    D.scale.values = zeros(D.Nchannels, D.Nevents);

    spm_progress_bar('Init', D.Nevents, 'Events downsampled'); drawnow;
    if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
    else, Ibar = [1:D.Nevents]; end

    for i = 1:D.Nevents
        d = squeeze(D.data(:, :, i));
        d2 = resample(d', Radc_new, D.Radc)';
        D.scale.values(:, i) = spm_eeg_write(fpd, d2, 2, D.datatype);
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
    D.events.start = round(D.events.start./(D.Radc/Radc_new));
    D.events.stop = size(d2, 2) - D.events.start - 1;
    D.Nsamples = size(d2, 2);
else
    % continuous
    D.scale.dim = 1;
    D.scale.values = zeros(D.Nchannels, 1);

    % adjust the timing information
    D.events.time = round(D.events.time*Radc_new/D.Radc);
    
    spm_progress_bar('Init', D.Nchannels, 'Channels downsampled'); drawnow;
    if D.Nchannels > 100, Ibar = floor(linspace(1, D.Nchannels, 100));
    else, Ibar = [1:D.Nchannels]; end

    for i = 1:D.Nchannels
        d = squeeze(D.data(i, :, :));
        d2 = resample(d, Radc_new, D.Radc);
        if i == 1
            data = zeros(D.Nchannels, length(d2));
        end
        data(i, :) = d2;      

        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
   
    end
    D.scale.values = spm_eeg_write(fpd, data, 2, D.datatype);
    D.Nsamples = size(d2, 1);

    
end

spm_progress_bar('Clear');

fclose(fpd);

D.events.time = round(D.events.time./(D.Radc/Radc_new));
D.Radc = Radc_new;
D.data = [];
D.fname = ['d' D.fname];

if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
