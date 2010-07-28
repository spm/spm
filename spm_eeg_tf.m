function [Dtf, Dtph] = spm_eeg_tf(S)
% Compute instantaneous power and phase in peri-stimulus time and frequency
% FORMAT [Dtf, Dtph] = spm_eeg_tf(S)
%
% S                     - input structure
%
% fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%   
%   S.channels          - cell array of channel names. Can include generic
%                         wildcards: 'All', 'EEG', 'MEG' etc.
%
%   S.frequencies      - vector of frequencies of interest
%
%   S.timewin          - time window of interest in PST in ms. 
%
%   S.method           - name for the spectral estimation to use. This
%                        corresponds to the name of a plug-in function that comes
%                        after 'spm_eeg_specest_' prefix.
%   S.settings         - plug-in specific settings
%
%   S.phase            - also save phase dataset (1) or not (0)
%                        phase dataset cannot be computed for some
%                        spectral estimation methods
% Output:
% Dtf                   - M/EEG object with power (also written on disk)
% Dtph                  - M/EEG object with phase (also written on disk)
%__________________________________________________________________________
% This is a modular function for which plugins can be developed implementing
% specific spectral estimation methods. There are 3 basic plugins presently
% implemented and they can be used as templates for new plugins.
% The name of a plugin function should start with 'spm_eeg_specest_'
%
% morlet (spm_eeg_specest_morlet) - Morlet wavelet transform
%
% hilbert (spm_eeg_specest_hilbert) - filtering + Hilbert transform
%
% ft_mtmconvol (spm_eeg_specest_ft_mtmconvol) - Fieldtrip implementation
%                                           of multi-taper spectral
%                                           analysis
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_tf.m 4021 2010-07-28 12:43:16Z vladimir $

SVNrev = '$Rev: 4021 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Time-Frequency'); spm('Pointer','Watch');

if nargin == 0
    S = [];
end

%-Ensure backward compatibility
%--------------------------------------------------------------------------
S = spm_eeg_compatibility(S, mfilename);

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

if isequal(D.type, 'continuous')
    error('Time-frequency analysis can only be applied to epoched data');
end

%-Configure the analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'channels')
    S.channels = 'all';
end

chanind = D.selectchannels(S.channels);

if isempty(chanind)
    error('No channels selected.');
end

if ~isfield(S, 'frequencies') || isempty(S.frequencies)
    S.frequencies = 1:48;
end

if ~isfield(S, 'timewin')
    S.timewin = 1e3*[D.time(1) D.time(end)];
end

timeind = D.indsample(1e-3*min(S.timewin)):D.indsample(1e-3*max(S.timewin));

if ~isfield(S, 'method')
    S.method = 'morlet';
    S.settings = [];
end

if ~isfield(S, 'phase')
    S.phase = 0;
end

if isfield(S, 'settings')
    S1 = S.settings;
else
    S1 = [];
end

S1.frequencies = S.frequencies;

%-Run the analysis on all trials
%--------------------------------------------------------------------------
spm_progress_bar('Init', D.ntrials, 'trials done');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

for k = 1:D.ntrials
    trial = feval(['spm_eeg_specest_' S.method], S1, D(chanind, timeind, k), D.time(timeind));
    
    if k == 1
        
        if isfield(trial, 'fourier')
            outdata = trial.fourier;
        else
            outdata = trial.pow;
        end
        
        Nchannels = size(outdata, 1);        
        Nfrequencies = size(outdata, 2);
        Nsamples = size(outdata, 3);
        
        %-Generate output datasets
        %--------------------------------------------------------------------------
        Dtf = clone(D, ['tf_' D.fnamedat], [Nchannels Nfrequencies Nsamples D.ntrials]);
        Dtf = Dtf.frequencies(:, trial.freq);
        Dtf = timeonset(Dtf, trial.time(1));
        Dtf = fsample(Dtf, 1/diff(trial.time(1:2)));
        Dtf = transformtype(Dtf, 'TF');
        
        Dtf = chanlabels(Dtf, 1:Nchannels, D.chanlabels(chanind));
        Dtf = badchannels(Dtf, 1:Nchannels, D.badchannels(chanind));
        Dtf = chantype(Dtf, 1:Nchannels, D.chantype(chanind));
        Dtf = coor2D(Dtf, 1:Nchannels, coor2D(D,chanind));
        
        if S.phase && isfield(trial, 'fourier')
            Dtph = clone(Dtf, ['tph_' D.fnamedat]);
            Dtph = transformtype(Dtph, 'TFphase');
        else
            if ~isfield(trial, 'fourier')
                warning('Phase cannot be estimated with the requested method. Estimating power only.');
            end
            
            Dtph = [];
        end
        
    end
    
    if isfield(trial, 'fourier')
        Dtf(:, :, :, k) = trial.fourier.*conj(trial.fourier);
        
        if S.phase
            Dtph(:, :, :, k) = angle(trial.fourier);
        end
    elseif isfield(trial, 'pow')
        Dtf(:, :, :, k) = trial.pow;
    else
        error('The plug-in returned unexpected output');
    end
        
    if ismember(k, Ibar), spm_progress_bar('Set', k); end
end

spm_progress_bar('Clear');

%-Save new M/EEG dataset(s)
%--------------------------------------------------------------------------
Dtf = Dtf.history('spm_eeg_tf', S);
save(Dtf);
if ~isempty(Dtph)
    Dtph = Dtph.history('spm_eeg_tf', S);
    save(Dtph);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency: done'); spm('Pointer','Arrow');
