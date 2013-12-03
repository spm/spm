function D = spm_eeg_downsample(S)
% Downsample M/EEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S               - optional input struct
% (optional) fields of S:
%   S.D           - MEEG object or filename of M/EEG mat-file
%   S.fsample_new - new sampling rate, must be lower than the original one
%   S.prefix     - prefix for the output file (default - 'd')
%
% D               - MEEG object (also written on disk)
%__________________________________________________________________________
% 
% This function uses the Signal Processing toolbox from The MathWorks:
%               http://www.mathworks.com/products/signal/
% (function resample.m) if present and an homebrew version otherwise
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 5773 2013-12-03 18:05:14Z guillaume $

SVNrev = '$Rev: 5773 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG downsampling'); spm('Pointer','Watch');

%-Test for the presence of Signal Processing Matlab toolbox
%--------------------------------------------------------------------------
flag_tbx = license('checkout','signal_toolbox') && ~isempty(ver('signal'));
if ~flag_tbx
    disp(['warning: using homemade resampling routine ' ...
          'as signal toolbox is not available.']);
end

if ~isfield(S, 'prefix'),       S.prefix = 'd';           end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

% This is to handle non-integer sampling rates up to a reasonable precision
P = round(10*S.fsample_new);
Q = round(10*D.fsample);

%-First pass: Determine new D.nsamples
%==========================================================================
if flag_tbx % Signal Proc. Toolbox
    nsamples_new = ceil(nsamples(D)*P/Q);
else
    d             = double(squeeze(D(1, :, 1)));
    [d2,alpha]    = spm_resample(d,P/Q);
    fsample_new   = D.fsample*alpha;
    S.fsample_new = fsample_new;
    disp(['Resampling frequency is ',num2str(fsample_new), 'Hz'])
    nsamples_new  = size(d2, 2);
end

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
    
    % now downsample blocks of channels
    chncnt=1;
    for blk=1:blknum
        spm_progress_bar('Set','ylabel','reading...');
        % load old meeg object blockwise into workspace
        blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
        Dtemp=D(blkchan,:,1);
        chncnt=chncnt+blksz;
        
        spm_progress_bar('Set','ylabel','writing...');
        if flag_tbx % Signal Proc. Toolbox
            Dnew(blkchan,:) = resample(Dtemp', P, Q)';
        else
            Dnew(blkchan,:) = spm_resample(Dtemp,P/Q);
        end
        spm_progress_bar('Set', blkchan(end));
    end    
else
    %-Epoched
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.ntrials, 'Trials downsampled'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
    else Ibar = 1:D.ntrials; end
    for i = 1:D.ntrials  
        if flag_tbx % signal proc. toolbox
            Dnew(:, :, i) = resample(spm_squeeze(D(:, :, i), 3)', P, Q)';
        else
            Dnew(:, :, i) = spm_resample(spm_squeeze(D(:, :, i), 3),P/Q);
        end       
        
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
        
    end
end;

spm_progress_bar('Clear');

%-Save new downsampled M/EEG dataset
%--------------------------------------------------------------------------
Dnew = fsample(Dnew, S.fsample_new);
D    = Dnew;
D    = D.history('spm_eeg_downsample', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG downsampling: done'); spm('Pointer','Arrow');


%==========================================================================
function [Y,alpha] = spm_resample(X,alpha)
% Basic resample function (when no Signal Proc. Toolbox)
% FORMAT Y = spm_resample(X,alpha)
% IN:
%   - X: a nXm matrix of n time series
%   - alpha: the ration of input versus output sampling frequencies. If
%   alpha>1, rs(X,alpha) performs upsampling of the time series.
% OUT:
%   - Y: nX[alpha*m] matrix of resampled time series
%   - alpha: true alpha used (due to rational rounding)
% This function operates on rows of a signal matrix. This means it can be
% used on a block of channels.

% Jean Daunizeau

N0     = size(X,2);
N      = floor(N0*alpha);
alpha  = N/N0;
Y      = fftshift(fft(X,[],2),2);
sy     = size(Y,2);
middle = floor(sy./2)+1;
if alpha>1 % upsample
    N2 = floor((N-N0)./2);
    if N0/2 == floor(N0/2)
        Y(:,1) = []; % throw away non symmetric DFT coef
    end
    Y  = [zeros(size(Y,1),N2),Y,zeros(size(Y,1),N2)];
else % downsample
    N2 = floor(N./2);
    Y  = Y(:,middle-N2:middle+N2);
end
Y      = alpha*ifft(ifftshift(Y,2),[],2);
