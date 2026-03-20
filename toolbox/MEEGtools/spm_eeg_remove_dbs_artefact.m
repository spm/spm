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
    S.timewin = [0.7 6.2];
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

% (LLM) The cleanest approach that guarantees C1 continuity (matching both value and slope at the join)
% and reaches exactly zero with zero slope is a cubic Hermite extension. It requires no toolboxes:

% --- Parameters ---
extend_s = 0.005;          % length of the appended tail in seconds
n_slope  = 100;             % samples used to estimate the terminal slope
fs       = D.fsample;
% --- Estimate endpoint value and slope ---
y_end   = template(end);
slope   = (template(end) - template(end - n_slope)) / (n_slope / fs);   % derivative at join

% --- Build the extension via cubic Hermite ---
% Boundary conditions:
%   t=0 : value = y_end,  derivative = slope
%   t=T : value = 0,      derivative = 0
n_ext = round(extend_s * fs);
t_ext = (1:n_ext).' / fs;          % time from the join point
T     = extend_s;

% Hermite basis functions (normalized interval u in [0,1])
u  = t_ext / T;
h00 =  2*u.^3 - 3*u.^2 + 1;       % value  at u=0
h10 =    u.^3 - 2*u.^2 + u;       % slope  at u=0 (scaled by T)
h01 = -2*u.^3 + 3*u.^2;           % value  at u=1
h11 =    u.^3 -   u.^2;           % slope  at u=1 (scaled by T, set to 0)

tail = h00 * y_end + h10 * (slope * T);   % h01*0 + h11*0 vanish

% --- Append to template ---
template_extended = [template tail'];
% end LLM

t_end = t(1)+length(template_extended)-1;

spm_progress_bar('Init', D.ntrials, 'trials interpolated');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end


for k = 1:D.ntrials
    for c = 1:length(indchannels)
        ep = D(indchannels(c), t(1):t(2), k);
        scale   = ep/template;
        D(indchannels(c), t(1):t_end, k) =  D(indchannels(c), t(1):t_end, k)-scale*template_extended;
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
