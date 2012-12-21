function D = spm_eeg_filter(S)
% Filter M/EEG data
% FORMAT D = spm_eeg_filter(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%   S.band    - filterband [low|high|bandpass|stop]
%   S.freq    - cutoff frequency(-ies) [Hz]
%
%  Optional fields:
%   S.type    - filter type [default: 'butterworth']
%                 'butterworth': Butterworth IIR filter
%                 'fir':         FIR filter (using MATLAB fir1 function)
%   S.order   - filter order [default: 5 for Butterworth]
%   S.dir     - filter direction [default: 'twopass']
%                 'onepass'         forward filter only
%                 'onepass-reverse' reverse filter only, i.e. backward in time
%                 'twopass'         zero-phase forward and reverse filter
%   S.prefix  - prefix for the output file [default: 'f']
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_filter.m 5160 2012-12-21 16:58:38Z guillaume $

SVNrev = '$Rev: 5160 $';


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG filter'); spm('Pointer', 'Watch');

if ~isfield(S, 'type'),         S.type   = 'butterworth'; end
if ~isfield(S, 'dir'),          S.dir    = 'twopass';     end
if ~isfield(S, 'prefix'),       S.prefix = 'f';           end
if ~isfield(S, 'order')
    if strcmp(S.type, 'butterworth')
        S.order = 5;
    else
        S.order = [];
    end
end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

%-Check band
%--------------------------------------------------------------------------
switch lower(S.band)
    
    case {'low','high'}
        if numel(S.freq)~=1
            error('Cutoff frequency should be a single number');
        end
        
        if S.freq < 0 || S.freq > D.fsample/2
            error('Cutoff must be > 0 & < half sample rate');
        end
        
    case {'bandpass','stop'}
        if S.freq(1) < 0 || S.freq(2) > D.fsample/2 || S.freq(1) > S.freq(2)
            error('Incorrect frequency band specification');
        end
        
    otherwise
        error('Incorrect filter band.')
end

%-Filter
%==========================================================================

%-Generate new meeg object with new filenames
Dnew = clone(D, [S.prefix fname(D)]);

%-Determine channels for filtering
Fchannels = D.indchantype('Filtered');

Fs = D.fsample;

if strcmp(D.type, 'continuous')
    %-Continuous data
    %----------------------------------------------------------------------
    spm_progress_bar('Init', nchannels(D), 'Channels filtered');
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
        % load meeg object blockwise into workspace
        blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
        if isempty(blkchan), break, end
        Dtemp=D(blkchan,:,1);
        chncnt=chncnt+blksz;
        %loop through channels
        for j = 1:numel(blkchan)
            
            if ismember(blkchan(j), Fchannels)
                Dtemp(j, :) = spm_eeg_preproc_filter(S, Dtemp(j,:), Fs);
            end
            
            if any(Ibar == j), spm_progress_bar('Set', blkchan(j)); end
            
        end
        
        % write Dtemp to Dnew
        Dnew(blkchan,:,1)=Dtemp;
        clear Dtemp;
    end
    
else
    %-Single trial or epoched
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.ntrials, 'Trials filtered');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = 1:D.ntrials; end
    
    for i=1:D.ntrials
        Dtemp = D(Fchannels, :, i);
        Dtemp = spm_eeg_preproc_filter(S, Dtemp, Fs);
        Dnew(Fchannels, 1:Dnew.nsamples, i) = Dtemp;
        
        if any(Ibar == i), spm_progress_bar('Set', i); end
    end
end

spm_progress_bar('Clear');

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName',''); spm('Pointer', 'Arrow');


%==========================================================================
% function dat = spm_eeg_preproc_filter(S, dat, Fs)
%==========================================================================
function dat = spm_eeg_preproc_filter(S, dat, Fs)

Fp  = S.freq;

if isequal(S.type, 'fir')
    type = 'fir';
else
    type = 'but';
end

N   = S.order;
dir = S.dir;

switch S.band
    case 'low'
        dat = ft_preproc_lowpassfilter(dat,Fs,Fp,N,type,dir);
    case 'high'
        dat = ft_preproc_highpassfilter(dat,Fs,Fp,N,type,dir);
    case 'bandpass'
        dat = ft_preproc_bandpassfilter(dat, Fs, Fp, N, type, dir);
    case 'stop'
        dat = ft_preproc_bandstopfilter(dat,Fs,Fp,N,type,dir);
end
