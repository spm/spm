function D = spm_eeg_interpolate_artefact(S)
% 'Baseline Correction' for M/EEG data
% FORMAT D = spm_eeg_interpolate_artefact(S)
%
% S        - optional input struct
% (optional) fields of S:
%   S.D    - MEEG object or filename of M/EEG mat-file with epoched data
%   S.time - 2-element vector with start and end of baseline period [ms]
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
%
% Subtract average baseline from all M/EEG and EOG channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_interpolate_artefact.m 6965 2016-12-08 13:47:06Z vladimir $

SVNrev = '$Rev: 6965 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG artefact interpolation'); spm('Pointer','Watch');

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


%-Give an error for TF data
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    error('This function is for time domain data only');
end

%-Get input parameters
%--------------------------------------------------------------------------
try
    time   = S.time;
catch
    time   = spm_input('Start and stop of artefact [ms]', '+1', 'i', '', 2);
    S.time = time;
end


%-Converting to sec
%--------------------------------------------------------------------------
time = time/1000;

%-Baseline Correction
%--------------------------------------------------------------------------
t(1) = D.indsample(time(1));
t(2) = D.indsample(time(2));

if t(1)<2 || t(2)>=D.nsamples
    error('There should be at least one sample before and after the artefact');
end

if any(isnan(t))
    error('The artefact segment was not defined correctly.');
end

indchannels = D.indchantype({'MEEG', 'EOG'});

S1         = [];
S1.D       = D;
S1.outfile = ['i' D.fname];
S1.updatehistory = 0;
D          = spm_eeg_copy(S1);

spm_progress_bar('Init', D.ntrials, 'trials interpolated');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end


for k = 1: D.ntrials
    for c = 1:length(indchannels)
        tmp = D(indchannels(c), [t(1)-1 t(1) t(2) t(2)+1], k);
        D(indchannels(c), t(1):t(2), k) = spline([t(1)-1 t(1) t(2) t(2)+1], tmp, t(1):t(2));
    end
    
    if ismember(k, Ibar), spm_progress_bar('Set', k); end
end

spm_progress_bar('Clear');

%-Update history
%--------------------------------------------------------------------------
D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG artefact interpolation: done'); spm('Pointer','Arrow');
