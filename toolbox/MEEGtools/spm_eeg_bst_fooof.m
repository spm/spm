function Dout = spm_eeg_bst_fooof(S)
% Remove the aperiodic component from the spectrum using the FOOOF algorithm
% Donoghue et al. (2020). Nature Neuroscience, 23, 1655-1665.
%
% This uses the Brainstorm implementation by Luc Wilson
%
% FORMAT  D = spm_eeg_bst_fooof(S)
%
% S         - struct (optional)
% (optional) fields of S:
% S.D - meeg object, filename or a list of filenames of SPM EEG files
% S.freq_range        - frequency range for fitting
% S.peak_width_limits - how wide the peaks can be
% S.max_peaks         - maximal number of peaks
% S.min_peak_height   - minimal peak height
% S.aperiodic_mode    - shape of the aperiodic component fixed|knee%
% S.peak_threshold    - threshold for detecting a peak
% S.peak_type         - Shape of the peak fit best|gaussian|cauchy
% S.line_noise_freq   - Line noise frequency 50|60Hz
% S.line_noise_width  - range around line noise peaks to interpolate
% S.guess_weight      - Parameter to weigh initial estimates during
%                       optimization none|weak|strong
% S.proximity_threshold  -  threshold to remove the smallest of two peaks
%                           if too close
%
% Output:
% D         - MEEG data struct with FOOOF-corrected spectra
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','FOOOF spectral correction',0);

if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select([1 inf], '\.mat$', 'Select M/EEG mat file(s)');
    S.D = D;
end

D = spm_eeg_load(D);

if ~(isequal(D.transformtype, 'TF') && isequal(D.type, 'evoked') && size(D, 3)==1)
    error('The input should be an averaged spectrum');
end

if ~isfield(S, 'freq_range')
    S.freq_range = spm_input('Frequency range','+1', 'r', '[5 45]', 2)';
end

if ~isfield(S, 'peak_width_limits')
    S.peak_width_limits = spm_input('Peak width limits','+1', 'r', '[4 10]', 2)';
end

if ~isfield(S, 'max_peaks')
    S.max_peaks = spm_input('Max num of peaks', '+1', 'n', '10', 1);
end

if ~isfield(S, 'min_peak_height')
    S.min_peak_height = spm_input('Min peak height', '+1', 'r', '0', 1);
end

if ~isfield(S, 'aperiodic_mode')
    S.aperiodic_mode = spm_input('Aperiodic mode?','+1','fixed|knee', char('fixed', 'knee'));
end

if ~isfield(S, 'peak_threshold')
    S.peak_threshold = spm_input('Peak threshold', '+1', 'r', '2', 1);
end

if ~isfield(S, 'peak_type')
    S.peak_type = spm_input('Peak type?','+1','best|gaussian|cauchy', char('best', 'gaussian', 'cauchy'));
end

if ~isfield(S, 'power_line')
    S.power_line = spm_input('Line noise frequency', '+1', '50Hz|60Hz', [50, 60]);
end

if ~isfield(S, 'guess_weight')
    S.guess_weight = spm_input('Guess weight?','+1','none|weak|strong', char('none', 'weak', 'strong'));
end

if ~isfield(S, 'proximity_threshold')
    S.proximity_threshold = spm_input('Proximity threshold', '+1', 'r', '2', 1);
end

S.thresh_after = true;
S.verbose = true;
S.return_spectrum = true;

hOT = exist('fmincon');

freq = D.frequencies;

toplot = 1;%max(D.nchannels, D.ntrials)<=8;
if toplot
    Fgraph   = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
end
k = 1;
peak_table = {};
for c = 1:D.nchannels
    for i = 1:D.ntrials
        
        pow = squeeze(D(c, :, :, i));
                 
        [fs, fg] = process_fooof('FOOOF_matlab', shiftdim(pow, -1), freq, S, hOT);
        
        if i==1 && c==1
            Dout = clone(D, ['F' D.fname], [D.nchannels length(fs) 1 D.ntrials], 0, 1);
            Dout = Dout.frequencies(':', fs);
        end
        
        Dout(c, :, 1, i) = fg.power_spectrum-log10(fg.ap_fit);
                
        if toplot
            %subplot(D.nchannels, D.ntrials, k)
            figure;
            plot(fs, log10(fg.fooofed_spectrum), 'r', 'LineWidth', 2.5);
            hold on
            plot(fs, log10(fg.ap_fit), 'b--', 'LineWidth', 2.5);
            plot(fs, fg.power_spectrum, 'k', 'LineWidth', 2.5);
        end
        
        if ~isequal(fg.peak_params, zeros(1, 3))
           peak_table = [peak_table; {char(D.chanlabels(c)), char(D.conditions(i)), fg.peak_params}];
        end
        
        k = k+1;
    end
end

Dout.peak_table = peak_table;
save(Dout);
