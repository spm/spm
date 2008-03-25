function D = spm_eeg_filter(S)
% filter M/EEG data
% FORMAT D = spm_eeg_filter(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% D       - filename of EEG-data file or EEG data struct
% filter  - struct with the following fields:
%    type       - type of filter, currently only 'butterworth'
%    PHz        - cutoff [Hz]
%    parameter  - filter coefficients
% 
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_filter filters M/EEG data.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_filter.m 1243 2008-03-25 23:02:44Z stefan $

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
    D = spm_eeg_load(D);
catch    
    error(sprintf('Trouble reading file %s', D));
end
    
try
    filter.type = S.filter.type;
catch
    filter.type =...
        spm_input('filter type', '+1', 'b', 'butterworth');
end

if strcmpi(filter.type,'butterworth')
    try
        filter.order = S.filter.order;
    catch
        filter.order = 5;
    end
    try
        filter.para = S.filter.para;
    catch
        filter.para = [];
    end
end

try
    filter.band  = S.filter.band ;
catch
    filter.band = cell2mat(...
        spm_input('filterband', '+1', 'm',...
        'lowpass|highpass|bandpass|stopband',...
         {'low','high','bandpass','stop'}));
end


try
    filter.PHz = S.filter.PHz;
catch
    if strcmp(filter.band, 'low')
        str = 'Cutoff [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r');
            if PHz > 0 && PHz < D.fsample/2, break, end
            str = 'Cutoff must be > 0 & < half sample rate';
        end
    elseif strcmp(filter.band, 'high')
        str = 'Cutoff [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r');
            if PHz > 0 && PHz < D.fsample/2, break, end
            str = 'Cutoff must be > 0 & < half sample rate';
        end
    elseif strcmp(filter.band, 'bandpass')
        str = 'band [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r', [], 2);
            if PHz(1) > 0 && PHz(1) < D.fsample/2 && PHz(1) < PHz(2), break, end
            str = 'Cutoff 1 must be > 0 & < half sample rate and Cutoff 1 must be < Cutoff 2';
        end
    elseif strcmp(filter.band, 'stop')
        str = 'band [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r', [], 2);
            if PHz(1) > 0 && PHz(1) < D.fsample/2 && PHz(1) < PHz(2), break, end
            str = 'Cutoff 1 must be > 0 & < half sample rate and Cutoff 1 must be < Cutoff 2';
        end
    end
    filter.PHz = PHz;
end

if strcmpi(filter.type, 'butterworth')
    % butterworth coefficients
    if(isempty(filter.para))
    [B, A] = butter(filter.order, 2*filter.PHz/D.fsample, filter.band);
        filter.para{1} = B;
        filter.para{2} = A;
    elseif(length(filter.para)~=2)
    errordlg('Need two parameters for Butterworth filter');
    end
else
    errordlg('Currently unknown filter type: %s', filter.type);
end

spm('Pointer', 'Watch');drawnow;

% generate new meeg object with new filenames
Dnew = newdata(D, ['f' fnamedat(D)], [D.nchannels D.nsamples D.ntrials], D.dtype);

spm_progress_bar('Init', D.ntrials, 'Events filtered'); drawnow;
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials

    d = squeeze(D(:, :, i));

    for j = 1:D.nchannels
        if strcmpi(filter.type, 'butterworth')
            d(j,:) = filtfilt(filter.para{1}, filter.para{2}, d(j,:));
        end
    end

    % base line correction
    d = d - repmat(mean(d(:, 1:indsample(D,0)), 2), 1, D.nsamples);
    %     Dnew = putdata(Dnew, 1:Dnew.nchannels, 1:Dnew.nsamples, i, d);
    Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end


spm_progress_bar('Clear');

save(Dnew);

spm('Pointer', 'Arrow');
