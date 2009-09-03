function [Dtf, Dtf2] = spm_eeg_tf(S)
% Compute instantaneous power and phase in peri-stimulus time and frequency
% FORMAT [Dtf, Dtf2] = spm_eeg_tf(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.tf               - structure with (optional) fields:
%     S.tf.frequencies - vector of frequencies (Hz)
%     S.tf.rm_baseline - baseline removal (yes/no: 1/0)
%     S.tf.Sbaseline   - 2-element vector: start and stop of baseline
%                        (if rm_baseline yes)
%     S.tf.Mfactor     - Morlet wavelet factor
%     S.tf.channels    - vector of channel indices for which to compute TF
%     S.tf.phase       - compute phase (yes/no: 1/0)
%     S.tf.pow         - compute power/magnitude: 1/0 [default: power]
%     S.tf.collchans   - collapse channels (yes/no: 1/0). Will average
%                        power over channels after power estimation.
%                        THIS OPTION HAS BEEN TEMPORARILY SWITCHED OFF.
% 
% Dtf                  - MEEG object with power data (also written to disk)
% Dtf2                 - MEEG object with phase data (also written to disk)
%                        if S.tf.phase == 1, empty array if not
%__________________________________________________________________________
%
% spm_eeg_tf estimates instantaneous power and phase of data using the
% continuous Morlet wavelet transform.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_tf.m 3350 2009-09-03 13:19:20Z vladimir $

SVNrev = '$Rev: 3350 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Time-Frequency'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, Dtf = []; Dtf2 = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

%-Get parameters
%--------------------------------------------------------------------------
try
    tf.frequencies = S.tf.frequencies;
catch
    tf.frequencies = ...
        spm_input('Frequencies (Hz)', '+1', 'r', '', [1, inf]);
    S.tf.frequencies = tf.frequencies;
end

try
    tf.rm_baseline = S.tf.rm_baseline;
catch
    tf.rm_baseline = ...
        spm_input('Removal of baseline', '+1', 'y/n', [1,0], 2);
    S.tf.rm_baseline = tf.rm_baseline;
end

if tf.rm_baseline
    try
        tf.Sbaseline = S.tf.Sbaseline;
    catch
        tf.Sbaseline = ...
            spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
        S.tf.Sbaseline = tf.Sbaseline;
    end
end

try
    tf.Mfactor = S.tf.Mfactor;
catch
    tf.Mfactor = ...
        spm_input('Which Morlet wavelet factor?', '+1', 'r', '7', 1);
    S.tf.Mfactor = tf.Mfactor;
end

try
    tf.channels = S.tf.channels;
catch
    tf.channels = ...
        spm_input('Select channels', '+1', 'i', num2str([1:D.nchannels]));
    S.tf.channels = tf.channels;
end

try
    tf.phase = S.tf.phase;
catch
    tf.phase = ...
        spm_input('Compute phase?', '+1', 'y/n', [1,0], 2);
    S.tf.phase = tf.phase;
end

try
    tf.pow = S.tf.pow;
catch
    tf.pow = 1;
    S.tf.pow = tf.pow;
end

% if length(tf.channels) > 1
%     try
%         tf.collchans = S.tf.collchans;    % 1 = collapse across channels, 0 = do not
%     catch
%         tf.collchans = spm_input('Collapse channels?', '+1', 'y/n', [1,0], 2);
%         S.tf.collchans = tf.collchans;
%     end
% else
%     tf.collchans = 0;
% end

%-Generate Morlet wavelets
%--------------------------------------------------------------------------
M = spm_eeg_morlet(tf.Mfactor, 1000/D.fsample, tf.frequencies);

% if tf.collchans
%     Nchannels = 1;
% else
    Nchannels = length(tf.channels);
% end

Nfrequencies = length(tf.frequencies);

Dtf = clone(D, ['tf1_' D.fnamedat], [Nchannels Nfrequencies D.nsamples D.ntrials]);
Dtf = Dtf.frequencies(:, tf.frequencies);

% fix all channels
sD = struct(Dtf);

for i = 1:length(tf.channels)
    lbl = D.chanlabels(tf.channels(i));
    sD.channels(i).label = lbl{1};
    
    if D.badchannels(tf.channels(i))
        sD.channels(i).bad = 1;
    end
    ctype = D.chantype(tf.channels(i));
    sD.channels(i).type = ctype{1};
%     Dtf = coor2D(Dtf, tf.channels(i), coor2D(D, tf.channels(i)));
    % units?

end
Dtf = meeg(sD);
Dtf = coor2D(Dtf, [1:length(tf.channels)], coor2D(D,tf.channels));

if tf.phase == 1
    Dtf2 = clone(D, ['tf2_' D.fnamedat], [Nchannels Nfrequencies D.nsamples D.ntrials]);
    Dtf2 = Dtf2.frequencies(:, tf.frequencies);
    Dtf2 = transformtype(Dtf2, 'TFphase');
    
    % fix all channels
    sD = struct(Dtf2);

    for i = 1:length(tf.channels)
        lbl = D.chanlabels(tf.channels(i));
        sD.channels(i).label = lbl{1};
        if D.badchannels(tf.channels(i))
            sD.channels(i).bad = 1;
        end
        ctype = D.chantype(tf.channels(i));
        sD.channels(i).type = ctype{1};
        % units?
    end
    Dtf2 = meeg(sD);
    %     Dtf2 = coor2D(Dtf2,[1:length(tf.channels)], coor2D(D,tf.channels));
    Dtf2 = coor2D(Dtf2, [1:length(tf.channels)], coor2D(D,tf.channels));
else
    Dtf2 = [];
end


spm_progress_bar('Init', D.ntrials, 'trials done');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = [1:D.ntrials]; end

for k = 1:D.ntrials

    d = zeros(length(tf.channels), Nfrequencies, D.nsamples);

    if tf.phase
        d2 = zeros(length(tf.channels), Nfrequencies, D.nsamples);
    end

    indat = squeeze(D(:, :, k));
    for j = 1:length(tf.channels)
        for i = 1 : Nfrequencies
            tmp = conv(indat(tf.channels(j), :), M{i});

            % time shift to remove delay
            tmp = tmp([1:D.nsamples] + (length(M{i})-1)/2);

            % power
            if tf.pow
                d(j, i, :) = tmp.*conj(tmp);
            else
                d(j, i, :) = abs(tmp);
            end

            % phase
            if tf.phase
                d2(j, i, :) = angle(tmp);
            end

        end
    end

%     if tf.collchans
%         d = mean(d, 1);
% 
%         if tf.phase
%             if tf.circularise
%                 tmp = cos(d2) + sqrt(-1)*sin(d2);
%                 d2 = double(abs(mean(tmp,1))./mean(abs(tmp),1));
%             else
%                 d2 = mean(d2,1);
%             end
%         end
%     end

    Dtf(:, :, :, k) = d;

    if tf.phase
        Dtf2(:, :, :, k) = d2;
    end

    if ismember(k, Ibar), spm_progress_bar('Set', k); end

end

spm_progress_bar('Clear');


%-Save new M/EEG dataset
%--------------------------------------------------------------------------
Dtf = Dtf.history('spm_eeg_tf', S);
save(Dtf);
if tf.phase
    Dtf2 = Dtf2.history('spm_eeg_tf', S);
    save(Dtf2);
end

%-Remove baseline over frequencies and trials
%--------------------------------------------------------------------------
if tf.rm_baseline == 1
    Dtf = spm_eeg_bc(struct('D',    Dtf, ...
                            'time', tf.Sbaseline,...
                            'save', false));
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency: done'); spm('Pointer','Arrow');
