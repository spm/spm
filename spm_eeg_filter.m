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
% $Id: spm_eeg_filter.m 865 2007-07-25 12:37:25Z christophe $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG filter setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

if ischar(D)
    P = spm_str_manip(D, 'H');
else
    P = D.path;
end


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

if strcmpi(D.filter.type,'butterworth')
    try
        D.filter.order = S.filter.order;
    catch
        D.filter.order = 5;
    end
    try
        D.filter.para = S.filter.para;
    catch
        D.filter.para = [];
    end
end

try
    D.filter.band  = S.filter.band ;
catch
    D.filter.band = cell2mat(...
        spm_input('filterband', '+1', 'm',...
        'lowpass|highpass|bandpass|stopband',...
         {'low','high','bandpass','stop'}));
end


try
    D.filter.PHz = S.filter.PHz;
catch
    if strcmp(D.filter.band, 'low')
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
    elseif strcmp(D.filter.band, 'high')
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
    elseif strcmp(D.filter.band, 'bandpass')
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
    elseif strcmp(D.filter.band, 'stop')
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
    if(isempty(D.filter.para))
	[B, A] = butter(D.filter.order, 2*D.filter.PHz/D.Radc, D.filter.band);
    	D.filter.para{1} = B;
    	D.filter.para{2} = A;
    elseif(length(D.filter.para)~=2)
	errordlg('Need two parameters for Butterworth filter');
    end
else
    errordlg('Currently unknown filter type: %s', D.filter.type);
end

spm('Pointer', 'Watch');drawnow;


% treat continuous and epoched data differently because of different
% scaling method

if size(D.data, 3) > 1
    % epoched

    % Prepare for writing data
    D.fnamedat = ['f' D.fnamedat];
    fpd = fopen(fullfile(P, D.fnamedat), 'w');

    d = zeros(D.Nchannels, D.Nsamples);
    D.scale = zeros(D.Nchannels, 1, D.Nevents);
    
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

        D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
    
    fclose(fpd);

else
    % continuous

    % make copy of the original file
    if ~copyfile(fullfile(P, D.fnamedat), fullfile(P, ['f' D.fnamedat]));
        error('Could not copy continuous file to %s', fullfile(P, ['f' D.fnamedat]));
    end

    % save existing matfile
    D.data = [];
    D.fnamedat = ['f' D.fnamedat];
    save(fullfile(P, ['f' D.fname]), 'D');

    % reload copied file
    D = spm_eeg_ldata(fullfile(P, ['f' D.fname]));
    

    D.scale = zeros(D.Nchannels, 1);

    spm_progress_bar('Init', D.Nchannels, 'Channels filtered'); drawnow;
    if D.Nchannels > 100, Ibar = floor(linspace(1, D.Nchannels, 100));
    else, Ibar = [1:D.Nchannels]; end

    for i = 1:D.Nchannels
        if strcmpi(D.filter.type, 'butterworth')
            d = filtfilt(D.filter.para{1}, D.filter.para{2}, squeeze(D.data(i, :, 1)));
        end

        D.scale(i) = spm_eeg_write(-1, d, 2, D.datatype);
        D.data(i, :) = d;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end

end

spm_progress_bar('Clear');


D.data = [];

D.fname = ['f' D.fname];
if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
