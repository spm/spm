function D = spm_eeg_downsample(S)
% function used for down-sampling EEG/MEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file
% fsample_new  - new sampling rate
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%_______________________________________________________________________
% This downsampling function uses 'resample' from the signal processing
% toolbox. Specify a new sampling rate, which must be lower than the
% original sampling rate. 
%_______________________________________________________________________
% Stefan Kiebel
% $Id: spm_eeg_downsample.m 2459 2008-11-12 11:27:23Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG downsample setup',0);

try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select EEG mat file');
    S.D = D;
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

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

spm('Pointer', 'Watch');drawnow;

% two passes

% 1st: Determine new D.nsamples
d = double(squeeze(D(1, :, 1)));
d2 = resample(d', P, Q)';
nsamples_new = size(d2, 2);

% generate new meeg object with new filenames
Dnew = clone(D, ['d' fnamedat(D)], [D.nchannels nsamples_new D.ntrials]);

% 2nd: resample all
if ~strcmp(D.type, 'continuous')
    spm_progress_bar('Init', D.ntrials, 'Events downsampled'); drawnow;
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
else
    spm_progress_bar('Init', D.nchannels, 'Channels downsampled'); drawnow;
    if D.nchannels > 100, Ibar = floor(linspace(1, D.nchannels,100));
    else Ibar = [1:D.nchannels]; end
end

now = clock;

for i = 1:D.ntrials
    for j = 1:D.nchannels
        d = double(squeeze(D(j, :, i)));
        d2 = resample(d', P, Q)';

        Dnew(j, 1:nsamples_new, i) = d2;

        if strcmp(D.type, 'continuous')
            if ismember(j, Ibar)
                spm_progress_bar('Set', j); drawnow;
            end
        end

    end
    if ~strcmp(D.type, 'continuous')
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end

end

elapsedTime = etime(clock,now);
fprintf(['Downsampling took ',num2str(elapsedTime),' sec\n'])

spm_progress_bar('Clear');

Dnew = putfsample(Dnew, fsample_new);
Dnew = putnsamples(Dnew, nsamples_new);

D = Dnew;
D = D.history('spm_eeg_downsample', S);
save(D);

spm('Pointer', 'Arrow');
