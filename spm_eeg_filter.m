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
% $Id: spm_eeg_filter.m 2389 2008-10-23 11:15:23Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG filter setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    S.D = D;
end

if ischar(D)
    P = spm_str_manip(D, 'H');
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
else
    P = D.path;
end


try
    filter.type = S.filter.type;
catch
    filter.type =...
        spm_input('filter type', '+1', 'b', 'butterworth');
    S.filter.type = filter.type;
end

if strcmpi(filter.type,'butterworth')
    try
        filter.order = S.filter.order;
    catch
        filter.order = 5;
        S.filter.order = filter.order;
    end
    try
        filter.para = S.filter.para;
    catch
        filter.para = [];
        S.filter.para = filter.para;
    end
end

try
    filter.band  = S.filter.band;
catch
    filter.band = cell2mat(...
        spm_input('filterband', '+1', 'm',...
        'lowpass|highpass|bandpass|stopband',...
         {'low','high','bandpass','stop'}));
     S.filter.band = filter.band;
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
            if PHz > 0 & PHz < D.fsample/2, break, end
            str = 'Cutoff must be > 0 & < half sample rate';
        end
    elseif strcmp(filter.band, 'high')
        str = 'Cutoff [Hz]';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [PHz, YPos] = spm_input(str, YPos, 'r', [], 1);
            if PHz > 0 & PHz < D.fsample/2, break, end
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
            if PHz(1) > 0 & PHz(1) < D.fsample/2 & PHz(1) < PHz(2), break, end
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
            if PHz(1) > 0 & PHz(1) < D.fsample/2 & PHz(1) < PHz(2), break, end
            str = 'Cutoff 1 must be > 0 & < half sample rate and Cutoff 1 must be < Cutoff 2';
        end
    end
    filter.PHz = PHz;
    S.filter.PHz = filter.PHz;
end

if strcmpi(filter.type, 'butterworth')
    % butterworth coefficients
    if(isempty(filter.para))
    [B, A] = butter(filter.order, filter.PHz/(D.fsample/2), filter.band);
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
Dnew = clone(D, ['f' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

% determine channels for filtering
Fchannels = unique([D.meegchannels, D.eogchannels, strmatch('LFP', D.chantype, 'exact')]);

ind = setdiff(1:D.nchannels, Fchannels);

if strcmp(D.type, 'continuous')

    % copy channels not to be filtered
    if ~isempty(ind)
        for i = 1:length(ind)
            Dnew(ind(i),:,1)=D(ind(i),:,1);
        end
    end
    
    spm_progress_bar('Init', length(Fchannels), 'Channels filtered'); drawnow;
    if length(Fchannels) > 100, Ibar = floor(linspace(1, length(Fchannels),100));
    else Ibar = [1:length(Fchannels)]; end

    for j = 1:length(Fchannels)
        
        if strcmpi(filter.type, 'butterworth')
            Dnew(Fchannels(j), :, 1) = filtfilt(filter.para{1}, filter.para{2}, squeeze(D(Fchannels(j), :, 1)));
        end

        if ismember(j, Ibar)
            spm_progress_bar('Set', j); drawnow;
        end

    end

else
    
    % copy channels not to be filtered
    if ~isempty(ind)
        Dnew(ind, :, 1:D.ntrials) = D(ind,:,1:D.ntrials);
    end   
    
    % single trial or epoched
    spm_progress_bar('Init', D.ntrials, 'Events filtered'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end


    for i = 1:D.ntrials

        d = squeeze(D(Fchannels, :, i));

        for j = 1:length(Fchannels)
            if strcmpi(filter.type, 'butterworth')
                d(j,:) = filtfilt(filter.para{1}, filter.para{2}, double(d(j,:)));
            end
        end

        % base line correction
        d = d - repmat(mean(d(:, 1:indsample(D,0)), 2), 1, D.nsamples);

        Dnew(Fchannels, 1:Dnew.nsamples, i) = d;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
end

spm_progress_bar('Clear');

D = Dnew;

% history
D = D.history('spm_eeg_filter', S);
save(D);

spm('Pointer', 'Arrow');

