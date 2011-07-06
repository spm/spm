function D = spm_eeg_downsample(S)
% Downsample M/EEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S               - optional input struct
% (optional) fields of S:
%   S.D           - MEEG object or filename of M/EEG mat-file
%   S.fsample_new - new sampling rate, must be lower than the original one
%   S.prefix      - prefix of generated file
%
% D               - MEEG object (also written on disk)
%__________________________________________________________________________
% 
% This function uses MATLAB Signal Processing Toolbox:
%               http://www.mathworks.com/products/signal/
% (function resample.m) if present and an homebrew version otherwise
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 4383 2011-07-06 15:55:39Z guillaume $

SVNrev = '$Rev: 4383 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG downsampling'); spm('Pointer','Watch');

%-Test for the presence of MATLAB Signal Processing Toolbox
%--------------------------------------------------------------------------
flag_TBX = license('checkout','signal_toolbox') & ~isdeployed;
if ~flag_TBX
    disp(['warning: using homemade resampling routine ' ...
          'as Signal Processing Toolbox is not available.']);
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
    fsample_new = S.fsample_new;
catch
    str = 'New sampling rate';
    YPos = -1;
    while 1
        if YPos == -1
            YPos = '+1';
        end
        [fsample_new, YPos] = spm_input(str, YPos, 'r');
        if fsample_new < D.fsample, break, end
        str = sprintf('Sampling rate must be less than original (%d)', round(D.fsample));
    end
    S.fsample_new = fsample_new;
end

try
    prefix = S.prefix;
catch
    prefix = 'd';
end

% This is to handle non-integer sampling rates up to a reasonable precision
P = round(10*fsample_new);
Q = round(10*D.fsample);

%-First pass: Determine new D.nsamples
%==========================================================================
if flag_TBX % Signal Proc. Toolbox
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
Dnew = clone(D, [prefix fnamedat(D)], [D.nchannels nsamples_new D.ntrials]);
t0 = clock;

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
    chncnt=1;
    for blk=1:blknum
        % load old meeg object blockwise into workspace 
        blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
        Dtemp=D(blkchan,:,1);
        chncnt=chncnt+blksz;
        
        %loop through channels
        for j = 1:numel(blkchan)
            d = Dtemp(j,1:D.nsamples);
            Dtemp(j,:)=0; % overwrite Dtemp to save memory
            if flag_TBX % Signal Proc. Toolbox
                Dtemp(j,1:nsamples_new) = resample(d', P, Q)';
            else
                Dtemp(j,1:nsamples_new) = spm_resample(d,P/Q);
            end
            spm_progress_bar('Set', blkchan(j));
        end

        % write Dtempnew to Dnew
        Dnew(blkchan,:,1)=Dtemp(:,1:nsamples_new,1);
        clear Dtemp

    end
    
else
    %-Epoched
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.ntrials, 'Trials downsampled'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        for j = 1:D.nchannels
            d = double(squeeze(D(j, :, i)));
            if flag_TBX % Signal Proc. Toolbox
                d2 = resample(d', P, Q)';
            else
                d2 = spm_resample(d,P/Q);
            end
            Dnew(j, 1:nsamples_new, i) = d2;
        end
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
    end
end

spm_progress_bar('Clear');

%-Display statistics
%--------------------------------------------------------------------------
fprintf('Elapsed time is %f seconds.\n',etime(clock,t0));               %-#

%-Save new downsampled M/EEG dataset
%--------------------------------------------------------------------------
Dnew = fsample(Dnew, (P/Q)*D.fsample);
D    = Dnew;
D    = D.history('spm_eeg_downsample', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG downsampling: done'); spm('Pointer','Arrow');


%==========================================================================
function [Y,alpha] = spm_resample(X,alpha)
% [Jean:] Basic resample function (when no Signal Proc. Toolbox)
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
