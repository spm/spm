function D = spm_eeg_remove_dbs_artefact(S)
% Template-based DBS artefact removal for EEG
% FORMAT D = spm_eeg_remove_dbs_artefact(S)
%
% S        - optional input struct
% (optional) fields of S:
%   S.D    - MEEG object or filename of M/EEG mat-file with epoched data
%   S.Da   - MEEG object or filename of M/EEG mat-file with artefact template
%   S.timewin - 2-element vector with start and end of baseline period [ms]
%   S.artchan - name of the channel containing artefact reference
%   S.artcond - condition name to take artefact from (if more than one)
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
% Based on Peeters et al.
% Electrophysiologic Evidence That Directional Deep Brain Stimulation Activates
% Distinct Neural Circuits in Patients With Parkinson Disease
% https://doi.org/10.1016/j.neurom.2021.11.002
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2026 Department of Imaging Neuroscience

if ~isfield(S, 'timewin')
    S.timewin = [0.7 3.5];
end


if ~isfield(S, 'artchan')
    S.artchan = 'Art';
end



%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename);
spm('FigName','DBS artefact removal'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

if ~isfield(S, 'Da')
    S.Da = S.D;
end

D = spm_eeg_load(D);


Da = spm_eeg_load(S.Da);

if ~isfield(S, 'artcond') && (Da.ntrials == 1)
    S.artcond = char(Da.condlist);
end

%-Give an error for TF data
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    error('This function is for time domain data only');
end

if ~isequal(D.type, 'evoked')
    error('This function works on evoked response');
end

time = S.timewin;

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

indchannels = D.indchantype('EEG');


if ~isequal(D.time, Da.time)
    error('Time axes mismatch');
end

S1         = [];
S1.D       = D;
S1.outfile = ['t' D.fname];
S1.updatehistory = 0;
D          = spm_eeg_copy(S1);

template = Da(Da.indchannel(S.artchan), t(1):t(2), Da.indtrial(S.artcond));


template_extended = append_forced_zero_tail(template, D.fsample, 0.15, 0.03)';

t_end = t(1)+length(template_extended)-1;

spm_progress_bar('Init', D.ntrials, 'trials interpolated');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end


for k = 1:D.ntrials
    for c = 1:length(indchannels)
        ep = D(indchannels(c), t(1):t(2), k);
        scale   = ep/template;
        D(indchannels(c), t(1):t_end, k) = D(indchannels(c), t(1):t_end, k)-scale*template_extended;                       
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

end

function [signal_ext, tail, info] = append_forced_zero_tail(signal, fs, extend_s, decayFracAtQuarter)
% Fast decay tail that is guaranteed to be zero at the end.
%
% Inputs:
%   signal             - Nx1 or 1xN waveform
%   fs                 - sampling rate (Hz)
%   extend_s           - extension length in seconds (e.g., 0.05 to 0.15)
%   decayFracAtQuarter - target remaining fraction at 25% of tail
%                        (e.g., 0.05 means 95% decayed by quarter-tail)
%
% Output:
%   signal_ext         - original + appended tail
%   tail               - appended segment
%   info               - parameters used

signal = signal(:);
N = numel(signal);

if nargin < 3 || isempty(extend_s)
    extend_s = 0.08; % 80 ms default
end
if nargin < 4 || isempty(decayFracAtQuarter)
    decayFracAtQuarter = 0.05; % aggressive/faster decay
end

nTail = max(2, round(extend_s * fs));
T = nTail / fs;
t = (1:nTail).' / fs;

% Start from actual endpoint (handles slight overshoot, positive or negative)
y0 = signal(end);

% Choose tau so raw exponential reaches decayFracAtQuarter at T/4
p = 0.25;
tau = -(p * T) / log(max(decayFracAtQuarter, 1e-6));

% Truncated/renormalized exponential:
% g(0)=1 and g(T)=0 exactly
e = exp(-t / tau);
eT = exp(-T / tau);
g = (e - eT) / (1 - eT);

tail = y0 * g;

% Guarantee exact zero at end sample
tail(end) = 0;

signal_ext = [signal; tail];

info = struct( ...
    'y0', y0, ...
    'extend_ms', 1000*T, ...
    'tau_ms', 1000*tau, ...
    'decayFracAtQuarter', decayFracAtQuarter, ...
    'nTail', nTail);
end