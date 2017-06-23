function D = spm_eeg_downsample(S)
% Downsample M/EEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S               - optional input struct
% (optional) fields of S:
%   S.D           - MEEG object or filename of M/EEG mat-file
%   S.method      - resampling method. Can be  'resample' [default],
%                   'decimate', 'downsample', 'fft'
%   S.fsample_new - new sampling rate, must be lower than the original one
%   S.prefix      - prefix for the output file [default: 'd']
%
% D               - MEEG object (also written on disk)
%__________________________________________________________________________
%
% This function uses the Signal Processing toolbox from The MathWorks:
%               http://www.mathworks.com/products/signal/
% (function resample.m) if present and spm_timeseries_resample.m otherwise.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG downsampling'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix = 'd';           end
if ~isfield(S, 'method'),       S.method = 'resample';    end

flag_tbx = license('checkout','signal_toolbox') && ~isempty(ver('signal'));
if ~flag_tbx && ~isequal(S.method, 'fft')
    S.method = 'fft';
    disp(['warning: switching to ''fft'' method ' ...
        'as signal toolbox is not available.']);
end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

% This is to handle non-integer sampling rates up to a reasonable precision
P = round(10*S.fsample_new)/10;
Q = round(10*D.fsample)/10;

%-First pass: Determine new D.nsamples
%==========================================================================
t             = ft_preproc_resample(D.time, Q, P, S.method);
nsamples_new  = size(t, 2);
fsample_new   = (nsamples_new/D.nsamples)*D.fsample;

if abs(S.fsample_new - fsample_new)<=0.1
    fsample_new = S.fsample_new;
else
    fsample_new = round(10*fsample_new)/10; 
end

fprintf('%-40s: %30s\n','Sampling frequency',[num2str(D.fsample),'Hz']); %-#
fprintf('%-40s: %30s\n','Resampling frequency',[num2str(fsample_new),'Hz']); %-#

%-Generate new meeg object with new filenames
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels nsamples_new D.ntrials]);

%-Second pass: resample all
%==========================================================================
if strcmp(D.type, 'continuous')
    %-Continuous
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.nchannels, 'Channels downsampled');
    
    % work on blocks of channels
    % determine block size, dependent on memory
    memsz  = spm('Memory');
    datasz = nchannels(D)*nsamples(D)*8; % datapoints x 8 bytes per double value
    blknum = ceil(datasz/memsz);
    blksz  = ceil(nchannels(D)/blknum);
    blknum = ceil(nchannels(D)/blksz);
    
    % now downsample blocks of channels
    chncnt = 1;
    for blk=1:blknum
        spm_progress_bar('Set','ylabel','reading...');
        % load old meeg object blockwise into workspace
        blkchan = chncnt:(min(nchannels(D), chncnt+blksz-1));
        Dtemp = D(blkchan,:,1);
        chncnt = chncnt+blksz;
        
        spm_progress_bar('Set','ylabel','writing...');
        
        Dnew(blkchan,:) = ft_preproc_resample(Dtemp, Q, P, S.method);
        
        spm_progress_bar('Set', blkchan(end));
    end
else
    %-Epoched
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.ntrials, 'Trials downsampled'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
    else Ibar = 1:D.ntrials; end
    for i = 1:D.ntrials
        Dnew(:, :, i) = ft_preproc_resample(spm_squeeze(D(:, :, i), 3), Q, P, S.method);
        
        if any(Ibar == i), spm_progress_bar('Set', i); end
        
    end
end

spm_progress_bar('Clear');

%-Save new downsampled M/EEG dataset
%--------------------------------------------------------------------------
Dnew = fsample(Dnew, S.fsample_new);
D    = Dnew;
D    = D.history('spm_eeg_downsample', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG downsampling: done'); spm('Pointer','Arrow');
