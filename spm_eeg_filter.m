function D = spm_eeg_filter(S)
% Filter M/EEG data
% FORMAT D = spm_eeg_filter(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%   S.filter  - struct with the following fields:
%      type       - type of filter, currently only 'butterworth'
%      band       - filterband [low|high|bandpass|stop]
%      PHz        - cutoff frequency [Hz]
%      parameter  - filter coefficients
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function filters M/EEG data and requires the signal processing 
% toolbox from The MathWorks:
%               http://www.mathworks.com/products/signal/
% (functions butter.m and filtfilt.m)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_filter.m 3401 2009-09-14 18:33:23Z guillaume $

SVNrev = '$Rev: 3401 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG filter'); spm('Pointer', 'Watch');

%-Test for the presence of required Matlab toolbox
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
    error('M/EEG filtering requires the Signal Processing Toolbox.');
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

%-Get parameters
%--------------------------------------------------------------------------
try
    filter.type = S.filter.type;
catch
    filter.type = spm_input('filter type', '+1', 'b', 'butterworth');
    S.filter.type = filter.type;
end

switch filter.type
    case 'butterworth'
        try
            filter.order   = S.filter.order;
        catch
            filter.order   = 5;
            S.filter.order = filter.order;
        end
        try
            filter.para    = S.filter.para;
        catch
            filter.para    = [];
            S.filter.para  = filter.para;
        end
    otherwise
        error('Unknown filter type.');
end

try
    filter.band   = S.filter.band;
catch
    filter.band   = cell2mat(...
        spm_input('filterband', '+1', 'm',...
        'lowpass|highpass|bandpass|stopband',...
        {'low','high','bandpass','stop'}));
    S.filter.band = filter.band;
end

try
    filter.PHz = S.filter.PHz;
catch
    switch lower(filter.band)
        case {'low','high'}
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
        case {'bandpass','stop'}
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
        otherwise
            error('unknown filter band.')
    end
    filter.PHz = PHz;
    S.filter.PHz = filter.PHz;
end

switch filter.type
    case 'butterworth'
        if isempty(filter.para)     
            [B, A] = butter(filter.order, filter.PHz/(D.fsample/2), filter.band);
            filter.para{1} = B;
            filter.para{2} = A;
        end
    otherwise
        error('Unknown filter type.');
end

%-
%--------------------------------------------------------------------------

% generate new meeg object with new filenames
Dnew = clone(D, ['f' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

% determine channels for filtering
Fchannels = unique([D.meegchannels, D.eogchannels]);

if strcmp(D.type, 'continuous')

    % continouous data
    spm_progress_bar('Init', nchannels(D), 'Channels filtered'); drawnow;
    if nchannels(D) > 100, Ibar = floor(linspace(1, nchannels(D),100));
    else Ibar = [1:nchannels(D)]; end

    % work on blocks of channels
    % determine blocksize
    % determine block size, dependent on memory
    memsz  = 2/3*spm('Memory');
    datasz = nchannels(D)*nsamples(D)*8; % datapoints x 8 bytes per double value
    blknum = ceil(datasz/memsz);
    blksz  = ceil(nchannels(D)/blknum);

    % now filter blocks of channels
    chncnt=1;
    for blk=1:blknum
        % load old meeg object blockwise into workspace
        blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
        if isempty(blkchan), break, end
        Dtemp=D(blkchan,:,1);
        chncnt=chncnt+blksz;
        %loop through channels
        for j = 1:numel(blkchan)

            if ismember(blkchan(j), Fchannels)
                Dtemp(j, :) = filtfilt(filter.para{1}, filter.para{2}, Dtemp(j,:));
            end

            if ismember(j, Ibar), spm_progress_bar('Set', blkchan(j)); end

        end

        % write Dtemp to Dnew
        Dnew(blkchan,:,1)=Dtemp;
        clear Dtemp;

    end;

else
    % single trial or epoched
    spm_progress_bar('Init', D.ntrials, 'Trials filtered'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end


    for i = 1:D.ntrials

        d = squeeze(D(:, :, i));

        for j = 1:nchannels(D)
            if ismember(j, Fchannels)
                d(j,:) = filtfilt(filter.para{1}, filter.para{2}, double(d(j,:)));
            end
        end

        Dnew(:, 1:Dnew.nsamples, i) = d;

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end
    
    disp('Baseline correction is no longer done automatically by spm_eeg_filter. Use spm_eeg_bc if necessary.');
end

spm_progress_bar('Clear');

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG filter: done'); spm('Pointer', 'Arrow');
