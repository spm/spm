function D = spm_eeg_bc(S)
% 'Baseline Correction' for M/EEG data
% FORMAT D = spm_eeg_bc(S)
%
% S        - optional input struct
% (optional) fields of S:
%   S.D    - MEEG object or filename of M/EEG mat-file with epoched data
%   S.time - 2-element vector with start and end of baseline period [ms]
%   S.save - save the baseline corrected data in a separate file [default: true]
%   S.updatehistory - update history information [default: true]
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
%
% Subtract average baseline from all M/EEG and EOG channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 3262 2009-07-09 12:10:53Z vladimir $

SVNrev = '$Rev: 3262 $';

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


%-Redirect to Time-Frequency baseline correction if necessary
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    D = spm_eeg_tf_rescale(S);
    return;
end

%-Get input parameters
%--------------------------------------------------------------------------
try
    time   = S.time;
catch
    time   = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
    S.time = time;
end


%-Converting to sec
%--------------------------------------------------------------------------
time = time/1000;

%-Baseline Correction
%--------------------------------------------------------------------------
t(1) = D.indsample(time(1));
t(2) = D.indsample(time(2));

if any(isnan(t))
    error('The baseline was not defined correctly.');
end

indchannels = [D.meegchannels D.eogchannels];


if ~isfield(S, 'save') || S.save
    S1         = [];
    S1.D       = D;
    S1.newname = ['b' D.fname];
    S1.updatehistory = 0;
    D          = spm_eeg_copy(S1);
end

spm_progress_bar('Init', D.ntrials, 'trials baseline-corrected');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end


for k = 1: D.ntrials
    tmp = mean(D(indchannels, t(1):t(2), k), 2);
    D(indchannels, :, k) = D(indchannels, :, k) - repmat(tmp, 1, D.nsamples);

    if ismember(k, Ibar), spm_progress_bar('Set', k); end
end

spm_progress_bar('Clear');

%-Save data
%--------------------------------------------------------------------------
if ~isfield(S, 'updatehistory') || S.updatehistory   
    D = D.history(mfilename, S);
    save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG baseline correction: done'); spm('Pointer','Arrow');
