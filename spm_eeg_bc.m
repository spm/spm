function D = spm_eeg_bc(S)
% 'Baseline Correction' for M/EEG data
% FORMAT D = spm_eeg_bc(S)
%
% S        - optional input struct
% (optional) fields of S:
%   S.D    - MEEG object or filename of M/EEG mat-file with epoched data
%   S.time - 2-element vector with start and end of baseline period [ms]
%   S.save - save the baseline corrected data on disk [default: false]
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
%
% Subtract average baseline from all M/EEG and EOG channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 2889 2009-03-17 12:02:04Z vladimir $

SVNrev = '$Rev: 2889 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG baseline correction'); spm('Pointer','Watch');

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

%-Get input parameters
%--------------------------------------------------------------------------
try
    time   = S.time;
catch
    time   = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
    S.time = time;
end

try 
    S.save;
catch
    S.save = false;
end

%-Converting to sec
%--------------------------------------------------------------------------
time = time/1000;

%-Baseline Correction
%--------------------------------------------------------------------------
t(1) = D.indsample(time(1));
t(2) = D.indsample(time(2));

indchannels = [D.meegchannels D.eogchannels];

spm_progress_bar('Init', D.ntrials, 'trials baseline-corrected');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

switch(transformtype(D))
    case {'TF','TFphase'}
        for k = 1: D.ntrials
            tmp = mean(D(:, :, t(1):t(2), k), 3);
            D(:, :, :, k) = D(:, :, :, k) - repmat(tmp, [1, 1, D.nsamples, 1]);

            if ismember(k, Ibar), spm_progress_bar('Set', k); end
        end

    case 'time'
        for k = 1: D.ntrials
            tmp = mean(D(indchannels, t(1):t(2), k), 2);
            D(indchannels, :, k) = D(indchannels, :, k) - repmat(tmp, 1, D.nsamples);
            
            if ismember(k, Ibar), spm_progress_bar('Set', k); end
        end

    otherwise
        error('Unknown transform type');
end

spm_progress_bar('Clear');

%-Save data
%--------------------------------------------------------------------------
if S.save
    error('Save option not handled yet.');
    %D = D.history('spm_eeg_bc', S);
    %save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG baseline correction: done'); spm('Pointer','Arrow');
