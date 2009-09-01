function D = spm_eeg_downsample(S)
% Downsample M/EEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S               - optional input struct
% (optional) fields of S:
%   S.D           - MEEG object or filename of M/EEG mat-file
%   S.fsample_new - new sampling rate, must be lower than the original one
%
% D               - MEEG object (also written on disk)
%__________________________________________________________________________
% 
% This function requires the Signal Processing toolbox from The MathWorks:
%               http://www.mathworks.com/products/signal/
% (function resample.m)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 3341 2009-09-01 14:23:49Z vladimir $

SVNrev = '$Rev: 3341 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG downsampling'); spm('Pointer','Watch');

%-Test for the presence of required Matlab toolbox
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
    error('M/EEG dowsampling requires the Signal Processing Toolbox.');
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

% This is to handle non-integer sampling rates up to a reasonable precision
P = round(10*fsample_new);
Q = round(10*D.fsample);

%-First pass: Determine new D.nsamples
%--------------------------------------------------------------------------
d = double(squeeze(D(1, :, 1)));
d2 = resample(d', P, Q)';
nsamples_new = size(d2, 2);

% generate new meeg object with new filenames
Dnew = clone(D, ['d' fnamedat(D)], [D.nchannels nsamples_new D.ntrials]);
t0 = clock;

%-Second pass: resample all
%--------------------------------------------------------------------------
if strcmp(D.type, 'continuous')
    spm_progress_bar('Init', D.nchannels, 'Channels downsampled');
    
    % work on blocks of channels
    % determine block size, dependent on memory
    try 
        evalc('memsz=2/3*feature(''memstats'');'); % 2/3 of largest block of contiguous memory, for Windows platforms
    catch
        memsz = 200*1024*1024; % 200 MB otherwise
    end
    datasz=nchannels(D)*nsamples(D)*8; % datapoints x 8 bytes per double value
    blknum=ceil(datasz/memsz);
    blksz=ceil(nchannels(D)/blknum);
    
    % now downsample blocks of channels
    chncnt=1;
    for blk=1:blknum
        % load old meeg object blockwise into workspace 
        blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
        Dtemp=D(blkchan,:,1);
        chncnt=chncnt+blksz;
        %loop through channels
        for j = 1:numel(blkchan)
            d = Dtemp(j,:);
            Dtemp(j,:)=0; % overwrite Dtemp to save memory
            Dtemp(j,1:nsamples_new) = resample(d', P, Q)';
            spm_progress_bar('Set', blkchan(j));
        end

        % write Dtempnew to Dnew
        Dnew(blkchan,:,1)=Dtemp(:,1:nsamples_new,1);
        clear Dtemp

    end
    
else
    spm_progress_bar('Init', D.ntrials, 'Trials downsampled'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        for j = 1:D.nchannels
            d = double(squeeze(D(j, :, i)));
            d2 = resample(d', P, Q)';

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
Dnew = fsample(Dnew, fsample_new);
D    = Dnew;
D    = D.history('spm_eeg_downsample', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG downsampling: done'); spm('Pointer','Arrow');
