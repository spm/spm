function D = spm_eeg_filter(S)
% Filter M/EEG data
% FORMAT D = spm_eeg_filter(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%   S.filter  - struct with the following fields:
%      type       - optional filter type, can be
%                    'but' Butterworth IIR filter (default)
%                    'fir' FIR filter using Matlab fir1 function
%      order      - filter order (default - 5 for Butterworth)
%      band       - filterband [low|high|bandpass|stop]
%      PHz        - cutoff frequency [Hz]
%      dir        - optional filter direction, can be
%                   'onepass'         forward filter only
%                   'onepass-reverse' reverse filter only, i.e. backward in time
%                   'twopass'         zero-phase forward and reverse filter
%                
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_filter.m 4127 2010-11-19 18:05:18Z christophe $

SVNrev = '$Rev: 4127 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG filter'); spm('Pointer', 'Watch');

if nargin == 0
    S = [];
end

%-Ensure backward compatibility
%--------------------------------------------------------------------------
S = spm_eeg_compatibility(S, mfilename);

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
if ~isfield(S, 'filter')
    S.filter = [];
end

if ~isfield(S.filter, 'band')
    S.filter.band   = cell2mat(...
        spm_input('filterband', '+1', 'm',...
        'lowpass|highpass|bandpass|stopband',...
        {'low','high','bandpass','stop'}));
end

if ~isfield(S.filter, 'type')
    S.filter.type = 'butterworth';
end

if ~isfield(S.filter, 'order')
    if strcmp(S.filter.type, 'butterworth')
        S.filter.order = 5;
    else
        S.filter.order = [];
    end
end

if ~isfield(S.filter, 'dir')
    S.filter.dir = 'twopass';
end


if ~isfield(S.filter, 'PHz')
    switch lower(S.filter.band)
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
    S.filter.PHz = PHz;
end

%-
%--------------------------------------------------------------------------

% generate new meeg object with new filenames
Dnew = clone(D, ['f' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

% determine channels for filtering
Fchannels = unique([D.meegchannels, D.eogchannels]);

Fs = D.fsample;

if strcmp(D.type, 'continuous')

    % continuous data
    spm_progress_bar('Init', nchannels(D), 'Channels filtered'); drawnow;
    if nchannels(D) > 100, Ibar = floor(linspace(1, nchannels(D),100));
    else Ibar = [1:nchannels(D)]; end

    % work on blocks of channels
    % determine blocksize
    % determine block size, dependent on memory
    memsz  = spm('Memory');
    datasz = nchannels(D)*nsamples(D)*8; % datapoints x 8 bytes per double value
    blknum = ceil(datasz/memsz);
    blksz  = ceil(nchannels(D)/blknum);
    blknum = ceil(nchannels(D)/blksz);

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
                Dtemp(j, :) = spm_eeg_preproc_filter(S.filter, Dtemp(j,:), Fs); 
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
                d(j,:) = spm_eeg_preproc_filter(S.filter,  double(d(j,:)), Fs);
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


%==========================================================================
function dat = spm_eeg_preproc_filter(filter, dat, Fs)

Fp  = filter.PHz;

if isequal(filter.type, 'fir')
    type = 'fir';
else
    type = 'but';
end

N   = filter.order;
dir = filter.dir;

switch filter.band
    case 'low'
        dat = ft_preproc_lowpassfilter(dat,Fs,Fp,N,type,dir);
    case 'high'
        dat = ft_preproc_highpassfilter(dat,Fs,Fp,N,type,dir);
    case 'bandpass'
        dat = ft_preproc_bandpassfilter(dat, Fs, Fp, N, type, dir);
    case 'stop'
        dat = ft_preproc_bandstopfilter(dat,Fs,Fp,N,type,dir);
end
