function D = spm_eeg_bc(S)
% 'Baseline Correction' for M/EEG data
% FORMAT D = spm_eeg_bc(S)
%
% S        - optional input struct
%      fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file with epoched data
%   S.timewin - 2-element vector with start and end of baseline period {ms}
%               [default: the negative times if present or the whole trial
%               otherwise]
%   S.save    - save the baseline corrected data in a separate file [default: true]
%   S.updatehistory - update history information [default: true]
%   S.prefix     - prefix for the output file [default: 'b']
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
%
% Subtract average baseline from all M/EEG and EOG channels
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 7132 2017-07-10 16:22:58Z guillaume $

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG baseline correction'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

%-Get input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'prefix'),          S.prefix = 'b';           end
if ~isfield(S, 'save'),            S.save = 1;               end
if ~isfield(S, 'updatehistory'),   S.updatehistory = 1;      end
if ~isfield(S, 'timewin')
    if D.time(1)<0
        timewin = [D.time(1) 0];
    else
        timewin = [D.time(1) D.time(end)];
    end
else
    timewin = 1e-3*S.timewin;
end

if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    error('Use spm_eeg_tf_rescale for TF data.');
end

%-Baseline Correction
%--------------------------------------------------------------------------
t(1) = D.indsample(timewin(1));
t(2) = D.indsample(timewin(2));

if any(isnan(t))
    error('The baseline was not defined correctly.');
end

indchannels = [D.indchantype('Filtered') D.indchantype('MEGCOMB')];

if S.save
    D = D.copy([S.prefix D.fname]);
end

spm_progress_bar('Init', D.ntrials, 'trials baseline-corrected');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end


for k = 1: D.ntrials
    tmp = mean(D(indchannels, t(1):t(2), k), 2);
    D(indchannels, :, k) = D(indchannels, :, k) - repmat(tmp, 1, D.nsamples);

    if any(Ibar == k), spm_progress_bar('Set', k); end
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
