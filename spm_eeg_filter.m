function D = spm_eeg_filter(S)
% low-pass filter EEG data
% FORMAT D = spm_eeg_filter(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% D		  - filename of EEG-data file or EEG data struct
% filter  - struct with the following fields:
%    type		- type of filter, currently only 'butterworth'
%    PHz        - cutoff [Hz]
%    parameter	- filter coefficients
% 
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_filter low-pass filters EEG/MEG epoched data.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_filter.m 149 2005-05-12 09:48:42Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG filter setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end
	
try
    D.filter.type = S.filter.type;
catch
	D.filter.type =...
	 	spm_input('filter type', '+1', 'b', 'butterworth');
end

try
    D.filter.band = S.filter.band;
catch
    D.filter.band =...
        spm_input('filterband', '+1', 'lowpass|bandpass');
end

try
    D.filter.PHz = S.filter.PHz;
catch
    if strcmp(D.filter.band, 'lowpass')
        str = 'Cutoff [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r');
            if PHz > 0 & PHz < D.Radc/2, break, end
            str = 'Cutoff must be > 0 & < half sample rate';
        end
    else
        % bandpass
        str = 'band [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r', [], 2);
            if PHz(1) > 0 & PHz(1) < D.Radc/2 & PHz(1) < PHz(2), break, end
            str = 'Cutoff 1 must be > 0 & < half sample rate and Cutoff 1 must be < Cutoff 2';
        end
    end
    D.filter.PHz = PHz;
end

if strcmpi(D.filter.type, 'butterworth')
    % butterworth coefficients
    D.filter.para = [];
    [B, A] = butter(5, 2*D.filter.PHz/D.Radc);
    D.filter.para{1} = B;
    D.filter.para{2} = A;
else
    error('Unknown filter type: %s', D.filter.type);
end

spm('Pointer', 'Watch');drawnow;

% Prepare for writing data
D.fnamedat = ['f' D.fnamedat];
fpd = fopen(fullfile(P, D.fnamedat), 'w');

% treat continuous and epoched data differently because of different
% scaling method

if size(D.data, 3) > 1
    % epoched
    d = zeros(D.Nchannels, D.Nsamples);
    D.scale.dim = [1 3];
    D.scale.values = zeros(D.Nchannels, D.Nevents);
    
    spm_progress_bar('Init', D.Nevents, 'Events filtered'); drawnow;
    if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
    else, Ibar = [1:D.Nevents]; end

    for i = 1:D.Nevents

        d = squeeze(D.data(:, :, i));

        for j = 1:D.Nchannels
            if strcmpi(D.filter.type, 'butterworth')
                d(j,:) = filtfilt(D.filter.para{1}, D.filter.para{2}, d(j,:));
            end
        end
        
        % base line correction
        d = d - repmat(mean(d(:,[1:abs(D.events.start)+1]),2), 1, D.Nsamples);

        D.scale.values(:, i) = spm_eeg_write(fpd, d, 2, D.datatype);
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
else
    % continuous
    D.scale.dim = 1;
    D.scale.values = zeros(D.Nchannels, 1);

    spm_progress_bar('Init', D.Nchannels, 'Channels filtered'); drawnow;
    if D.Nchannels > 100, Ibar = floor(linspace(1, D.Nchannels, 100));
    else, Ibar = [1:D.Nchannels]; end

    for i = 1:D.Nchannels
        if strcmpi(D.filter.type, 'butterworth')
            d = filtfilt(B, A, squeeze(D.data(i, :, :)));
        end
        
        if i == 1
            data = zeros(D.Nchannels, D.Nsamples);
        end
        data(i, :) = d;
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
    D.scale.values = spm_eeg_write(fpd, data, 2, D.datatype);

end

spm_progress_bar('Clear');

fclose(fpd);

D.data = [];

D.fname = ['f' D.fname];
if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
