function [cal_D, positions] = spm_opm_calibrate_from_coils(S)
%   [cal_D, positions] = spm_opm_calibrate_from_coils(S) estimates 
%   the 3D positions, orientations, and gains of magnetometers based on 
%   recordings during external coil activation. The function fits observed 
%   signals to a model of magnetic dipoles using both amplitude and phase 
%   errors. Results are then used to calibrate the magnetometer array.
%            
%   Input:
%       S : structure with fields
%           - D:               MEEG D object to be calibrated (required)
%           - freq:            Coil signal frequency in Hz (default: 10)
%           - chan_name_format Format of sensor/axis for channel (default:
%                              '^(X|Y|Z)(.*?)$'
%           - amplitude_range: 1x2 [min max] amplitude range (fT) for filtering
%                              (default: [0.1 1000]*1e3)
%           - coil_properties: structure with fields:
%               positions:     n x 3 dipole positions (m)
%               orientations:  n x 3 dipole unit vector
%               resistance:    default, 2120
%               inner_r:       default, 3.68/1000 (m)
%               outer_r:       default, 10.94/1000 (m)
%               layers:        default, 6
%               layer_turns:   default, 45
%           - t_delay:         timing offset for phase correction (s), 
%                              default: 'auto')
%           - phase_tol:       tolerance for phase filter applied (rads), 
%                              default: 'auto')
%           - coil_trigger:    structure with fields:
%               event_channel: channel name for triggers (default: 'A8')
%               event_design:  order of coil operations. (default: 
%                              repmat(1:16,[1 4]). Index with respect to
%                              coil_properies.
%               event_duration:in seconds (default: 1);
%               event_trim:    1 x 2, seconds to trim from [start, end].
%                              (default: [1/S.freq, 1/S.freq])
%               signal_channel:channel name for coil signals (default: 'A16')
%           - plot_output:     boolean
%           - highpass:        highpass filter freq (Hz), default 5
%           - channels:        cell array of channels. If empty, include all
%           - use_similarity   boolean, whether to use starting parameters
%                              based on similar observations across
%                              channels, default: true;
%           - fast             boolean, skip permutations where possible
%                              due to good starting conditions. Set to
%                              false if bad channels are missed, (default:
%                              true)
%           - balance          boolean, update existing forward model/s, 
%                              (default: false)
%           - estimation: structure with fields:
%               permutations:  number of permutations to run (default: 6)
%               min_ops:       minimum number of valid operations to
%                              estimate axis (default: 8)
%               sub_sel_num:   number of operations to estimate for initial 
%                              estimate/s (default: 5)
%               phase_tol:     tolerance for phase filter (default: 'auto')
%               pos_bounds:    lower/upper bounds per dimension (m) (default:
%                              [-0.15, 0.15; -0.15, 0.15; -0.3, -0.05]
%               gain_bounds:   [0.666, 1.5];
%               parallel:      boolean, use parallel processing, default
%                              false
%           - output_filename: filename for saving table;
%
%   Output:
%       cal_D:                 calibrated D object (spm_opm_apply_calibration)
%       positions: table containing calibrated sensor information
%           - name: sensor/channel name
%           - Px, Py, Pz: position coordinates (typically in metres)
%           - Ox, Oy, Oz: orientation vector components
%           - gain: channelwise gain factor
%           - cal_error: calibration error (RMS)
%
%   Example usage:
%       % Load file with HALO activation
%       S = [];
%       S.data = 'D:\data\sub-1_ses-1_task-halo_run-1_meg.lvm';
%       D = spm_opm_create(S);
%       % Calibrate
%       S = [];
%       S.D = D;
%       [cal_D, positions] = spm_opm_calibrate_from_coils(S);
%       % Load data from experiment (ensure no recalibration after HALO)
%       S = [];
%       S.data = 'D:\data\sub-1_ses-1_task-experiment_run-1_meg.lvm';
%       D2 = spm_opm_create(S);
%       % Apply calibration to data from experiment
%       S = [];
%       S.D = D2;
%       S.positions = positions;
%       [cal_D2] = spm_opm_apply_calibration(S)
%   
%   Notes:
%       The cal_error is reported per axis, but the optimisation is based on 
%           a combined cal_error term for all 'good' axes in a sensor.
%       cal_error is the RMS of amplitude error between modelled and observed 
%           responses. No threshold is set other than successful fit, but 
%           values greater than ~0.1 may be a reasonable threshold for a bad 
%           channel.
%       See Hill et al. (2025) https://doi.org/10.1162/imag_a_00535 for 
%           details on the HALO system (Quspin Inc. CO, USA) used as the 
%           default parameters in this function.
%       Several steps are taken to ensure only good data goes into the
%           estimation, including checking phase consistency, amplitude
%           correlation, amplitude thresholding, and removal of axes not
%           contributing to the fit.
%
% Author: Nicholas Alexander (n.alexander@ucl.ac.uk)
% Copyright: Department of Imaging Neuroscience, UCL, 2026
%
% TODO add precalculated library lookup of positions around sphere for phase
%% Setup
S = input_read_out(S);
amplitude_range = S.amplitude_range;
fs = S.D.fsample;

% Ideally these would come from the data/inputs
num_phase_bins = 1000;
min_valid_for_similar = 10;
good_similarity_threshold = 0.7;
corr_threshold = 0.9;
axis_plotting_colours = [252, 141, 98; 102, 194, 165; 141, 160, 203] ./ 255;

%% Coil coregistration
%  Coil information
coil_properties = S.coil_properties;
freq = S.freq;
t_delay = S.t_delay;
phase_tol = S.phase_tol;

% Moment values
turn_width = (coil_properties.outer_r - coil_properties.inner_r) / coil_properties.layer_turns;
turn_indices = 1:coil_properties.layer_turns;
turn_r = coil_properties.inner_r + ((turn_indices - 0.5) .* turn_width);
turn_area = repmat(pi .* (turn_r .^ 2), [1,coil_properties.layers]);
total_coil_area = sum(turn_area);

%% Highpass
if ~isempty(S.highpass)
    S_r = [];
    S_r.D = S.D;
    S_r.freq = S.highpass;
    S_r.band = 'high';
    S_r.prefix = 'proc_';
    proc_D = spm_eeg_ffilter(S_r);
    proc_D.save();
end

%% Process trigger information
% Trigger channel
trig_chan = S.coil_trigger.event_channel;
trig_idx  = S.D.indchannel(trig_chan);
trig_dat = squeeze(S.D(trig_idx,:,:));

% Coil operations
op_start = find(diff(trig_dat > 0.5) == 1) + 1;
if numel(op_start) ~= numel(S.coil_trigger.event_design)
    warning(['Number of trigger operations does not match input design. ' ...
        'Attempting to crop, but you should check this'])
    if numel(op_start) > numel(S.coil_trigger.event_design)
        op_start = op_start(1:numel(S.coil_trigger.event_design));
    elseif numel(op_start) < numel(S.coil_trigger.event_design)
        S.coil_trigger.event_design = S.coil_trigger.event_design(numel(op_start));
    end
end
op_end = op_start + round(S.coil_trigger.event_duration * fs) - 1;

% Optional trimming
if isfield(S.coil_trigger,'event_trim') 
    op_start  = op_start + round(S.coil_trigger.event_trim(1) * fs);
    op_end    = op_end - round(S.coil_trigger.event_trim(2) * fs);
    bounds = [op_start(1)-S.coil_trigger.event_trim(1)*fs, op_end(end)+S.coil_trigger.event_trim(2)*fs];
else
    bounds = [op_start(1), op_end(end)];
end

num_ops = length(S.coil_trigger.event_design);

% This channel contains the signal sent to the coils
coil_chan = S.coil_trigger.signal_channel;
coil_dat = S.D(S.D.indchannel(coil_chan),:,:);

%% Process channels
% Select channels
meg_idx = find(ismember(S.D.chantype,'MEGMAG'));
all_chans = S.D.chanlabels(meg_idx);

if ~isempty(S.channels)
    keep_mask = ismember(all_chans, S.channels);
    meg_idx = meg_idx(keep_mask);
    all_chans = all_chans(keep_mask);
end

% Axis info
in_parentheses = regexp(S.chan_name_format, '\(([^)]+)\)', 'tokens');
axis_info_idx = find(~cellfun(@isempty, regexp([in_parentheses{:}], '\|')));
axis_group_idx = axis_info_idx;
axes_labels = strsplit(in_parentheses{axis_info_idx}{1}, '|');
sensor_info = struct();

% Extract sensor names from channels
all_sensor_names = cell(1,numel(all_chans));
for chan_idx = 1:numel(all_chans)
    chan = all_chans{chan_idx};
    label_tokens = regexp(chan, S.chan_name_format, 'tokens', 'once');
    if isempty(label_tokens)
        continue
    end
    sensor_name_idx = setdiff(1:numel(label_tokens), axis_group_idx);
    sensor_name = strjoin(label_tokens(sensor_name_idx), '');
    all_sensor_names{chan_idx} = sensor_name;
end
sens_names = unique(all_sensor_names(~cellfun(@isempty, all_sensor_names)));
[~,sens_sort_idx] = sort(cellfun(@str2double, sens_names));
sens_names = sens_names(sens_sort_idx);

% Channel info per sensor
for sens_idx = 1:numel(sens_names)
    sensor_info(sens_idx).sensor_name = sens_names{sens_idx};
    for axis_idx = 1:numel(axes_labels)
        sensor_info(sens_idx).axis.(axes_labels{axis_idx}) = [];
    end
end

% Fill in channel indices per axis
for chan_idx = 1:numel(all_chans)
    chan = all_chans{chan_idx};
    label_tokens = regexp(chan, S.chan_name_format, 'tokens', 'once');
    if isempty(label_tokens)
        continue
    end

    axis_label = label_tokens{axis_group_idx};
    sensor_name_idx = setdiff(1:numel(label_tokens), axis_group_idx);
    sensor_name = strjoin(label_tokens(sensor_name_idx), '');

    % Find sensor index
    sens_idx = find(strcmp({sensor_info.sensor_name}, sensor_name), 1);
    if isempty(sens_idx)
        continue
    end

    % Store channel index in the appropriate axis field
    sensor_info(sens_idx).axis.(axis_label) = chan_idx;
end

num_axis = length(axes_labels);
num_sens = numel(sensor_info);

%% Plot example input data
sub_h = cell(1,2 + num_axis);
if S.plot_output
    % Visualise
    signal_figure_handle = figure;
    sub_h{1} = subplot(2 + num_axis,1,1);
    plot(S.D.time(bounds(1):bounds(2)), trig_dat(bounds(1):bounds(2)), 'Color', [231,138,195] ./ 255)
    y_lims = [min(trig_dat(bounds(1):bounds(2))), max(trig_dat(bounds(1):bounds(2)))];
    ylim(y_lims + [-range(y_lims), range(y_lims)] * 0.1)
    title('Coil timing')
    for op_idx = 1:num_ops
        x_pos = mean([op_start(op_idx), op_end(op_idx)]) ./ fs;
        y_pos = mean(trig_dat(op_start:op_end)) * 0.8;
        text(x_pos, y_pos, num2str(S.coil_trigger.event_design(op_idx)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontWeight','bold', 'FontSize',5, 'Color','k');
    end
   
    sub_h{2} = subplot(2 + num_axis,1,2);
    plot(S.D.time(bounds(1):bounds(2)), coil_dat(bounds(1):bounds(2)),'Color', [166,216,84]./255)
    y_lims = [min(coil_dat(bounds(1):bounds(2))), max(coil_dat(bounds(1):bounds(2)))];
    ylim(y_lims + [-range(y_lims),range(y_lims)]*0.1)
    title('Coil signal (V)')
end

%% Estimate sensor calibrations
% Plot coil array
if S.plot_output
    figure_handle = figure;
    hold on
    tmp_pos = coil_properties.positions;
    tmp_ori = coil_properties.orientations(1,:);
    num_coils = length(tmp_pos(:,1));
    for coil_idx = 1:num_coils
        quiver3(tmp_pos(coil_idx, 1), tmp_pos(coil_idx, 2), tmp_pos(coil_idx, 3), ...
            tmp_ori(1), tmp_ori(2), tmp_ori(3),0.05, 'k');
        text(tmp_pos(coil_idx, 1), tmp_pos(coil_idx, 2), tmp_pos(coil_idx, 3),num2str(coil_idx))
    end
end

% Precalculate phase and amplitude for coils
% Select all operations
op_coil_positions = zeros(num_ops, 3);
op_coil_orientations = zeros(num_ops, 3);
op_coil_amplitudes = zeros(num_ops,1);
op_coil_phases = zeros(num_ops,1);
for op_idx = 1:num_ops
    coil_idx = S.coil_trigger.event_design(op_idx);
    op_coil_orientations(op_idx, 1:3) = coil_properties.orientations(coil_idx,1:3);
    op_coil_positions(op_idx, 1:3) = coil_properties.positions(coil_idx,1:3);
    
    % Phase at specified frequency
    op_coil_signal = coil_dat(op_start(op_idx):op_end(op_idx));
    num_samples = length(op_coil_signal);
    time = 0:1/fs:((num_samples-1)/fs);
    complex_sum = sum(op_coil_signal .* exp(-1i * 2 * pi * freq * time), 2);
    op_coil_phases(op_idx) = angle(complex_sum);
    op_coil_amplitudes(op_idx) = 2 * abs(complex_sum) / num_samples;
end

% Dipole moment calculation
op_coil_currents = op_coil_amplitudes ./ coil_properties.resistance;
op_coil_moments = total_coil_area .* op_coil_currents;

% Repeat for sensors
observed_amplitudes = nan(num_sens, num_axis, num_ops);
observed_phases = nan(num_sens, num_axis, num_ops);
phase_diff = nan(num_sens, num_axis, num_ops);
for sens_idx = 1:num_sens
    % Select all channels from sensor
    cur_sensor_info = sensor_info(sens_idx);
    cur_chan_indices = nan(1,num_axis);
    for axis_idx = 1:num_axis
        axis_label = axes_labels{axis_idx};
        if isfield(cur_sensor_info.axis, axis_label) && ~isempty(cur_sensor_info.axis.(axis_label))
            cur_chan_indices(axis_idx) = cur_sensor_info.axis.(axis_label);
        end
    end

    valid_axis = find(~isnan(cur_chan_indices));

    if ~isempty(S.highpass)
        sens_dat = proc_D(meg_idx(cur_chan_indices(valid_axis)),:,:);
    else
        sens_dat = S.D(meg_idx(cur_chan_indices(valid_axis)),:,:);
    end
    
    % Prepare observed signals
    op_num_samples = length(op_start(1):op_end(1));
    observed_signals = nan(num_axis, num_ops, op_num_samples);
    for op_idx = 1:num_ops
        for valid_axis_idx = 1:length(valid_axis)
            observed_signals(valid_axis(valid_axis_idx), op_idx, 1:op_num_samples) = sens_dat(valid_axis_idx, op_start(op_idx):op_end(op_idx));
        end
    end

    % Calculate observed amplitudes
    for valid_axis_idx = 1:length(valid_axis)
        for op_idx = 1:num_ops
            % Select data
            observed_signal = squeeze(observed_signals(valid_axis(valid_axis_idx), op_idx,:));
    
            % Phase and amplitude at specified frequency
            num_samples = length(observed_signal);
            time = 0:1/fs:((num_samples-1)/fs);
            complex_sum = sum(observed_signal' .* exp(-1i * 2 * pi * freq * time), 2);
            observed_phases(sens_idx, valid_axis(valid_axis_idx), op_idx) = angle(complex_sum);
            observed_amplitudes(sens_idx, valid_axis(valid_axis_idx), op_idx) = 2 * abs(complex_sum) / num_samples;

            phase_diff(sens_idx, valid_axis(valid_axis_idx), op_idx) = angle(exp(1j * (observed_phases(sens_idx, valid_axis(valid_axis_idx), op_idx) - op_coil_phases(op_idx))));

        end
    end
end

% Phase delay/tolerance estimation
figure
hold on
if ischar(t_delay) && strcmp(t_delay, 'auto')
    % Wrapping
    phase_diff_all = [phase_diff; phase_diff + pi];
    phase_diff_all = phase_diff_all(:);
    phase_diff_all(phase_diff_all < -0.5*pi) = [];
    phase_diff_all(phase_diff_all > 1.5*pi) = [];
    phase_diff_all(isnan(phase_diff_all)) = [];

    % Identify most common delay (phase and anti-phase)
    [counts, edges] = histcounts(phase_diff_all, num_phase_bins);
    background_level = median(counts);
    noise_level = std(counts);
    prom_thresh = background_level + 3 * noise_level;
    [~, peak_x, peak_w] = findpeaks(counts, 'NPeaks', 2, 'MinPeakProminence', prom_thresh);
    
    % Get phase values
    bin_centres = edges(1:end-1) + diff(edges)/2;
    peak_phases = bin_centres(peak_x);
    
    if isempty(peak_phases)
        warning('Auto phase-delay detection failed')
        t_delay = 0;
    else
        if numel(peak_phases) == 1
            dominant_phase = peak_phases(1);
            t_delay = (dominant_phase / (2*pi*freq)) * 1000;
        elseif numel(peak_phases) >= 2
            peak_phase_diff = abs(peak_phases(2) - peak_phases(1));
            phase_antiphase_tol = 0.1;
            if abs(peak_phase_diff - pi) < phase_antiphase_tol
                dominant_phase = angle(exp(1j * mean([peak_phases(1), peak_phases(2) - pi])));
                t_delay = (dominant_phase / (2*pi*freq)) * 1000;
            else
                warning('Auto phase-delay detection failed')
                t_delay = 0;
            end
        end
    end
    
    fprintf('Detected phase/anti-phase pair with estimated delay: %.2f ms\n', t_delay);
    histogram(phase_diff, num_phase_bins,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[102,194,165]./255)

    % Phase tolerance
    if ischar(phase_tol) && strcmp(phase_tol, 'auto')
        bin_width = edges(2) - edges(1);
        phase_tol = max(peak_w * bin_width);
        fprintf('Setting phase tolerance to: %.2f rads\n', phase_tol);
    end
end

% Update phase difference
for sens_idx = 1:num_sens
    for axis_idx = 1:num_axis
        for op_idx = 1:num_ops
            observed_phases(sens_idx, axis_idx, op_idx) = observed_phases(sens_idx, axis_idx, op_idx) - 2*pi*freq*(t_delay/1000);
            phase_diff(sens_idx, axis_idx, op_idx) = angle(exp(1j * (observed_phases(sens_idx, axis_idx, op_idx) - op_coil_phases(op_idx))));
        end
    end
end
histogram(phase_diff, num_phase_bins,'DisplayStyle','stairs','LineWidth',1,'EdgeColor',[252,141,98]./255)
legend({'Original','Phase delay corrected'})
title(['Phase delay of ', num2str(t_delay),  ' ms applied'])
xlabel('Phase difference from coil to sensor')
ylabel('Density of observations')
hold off

% Phase relation filter
phase_relation = nan(size(observed_phases));
for sens_idx = 1:num_sens
    for axis_idx = 1:num_axis
        for op_idx = 1:num_ops
            cur_phase_diff = phase_diff(sens_idx, axis_idx, op_idx);
            if abs(cur_phase_diff) < phase_tol
                phase_relation(sens_idx, axis_idx, op_idx) = 1;
            elseif abs(abs(cur_phase_diff) - pi) < phase_tol
                phase_relation(sens_idx, axis_idx, op_idx) = -1;
            else
                phase_relation(sens_idx, axis_idx, op_idx) = NaN;
            end
        end
    end
end

% Amplitude correlation check
[~,~,position_group] = unique(op_coil_positions,'rows');
num_unique_pos = max(position_group);
grouped_amplitude_corr = nan(num_sens,num_axis);
for sens_idx = 1:num_sens
    for pos_idx = 1:num_unique_pos
        op_idx = find(position_group == pos_idx);
        for axis_idx = 1:num_axis
            cur_observed_amplitudes = squeeze(observed_amplitudes(sens_idx,axis_idx, op_idx));
            cur_coil_amplitudes = op_coil_moments(op_idx);
            current_amplitude_corr = corr(cur_observed_amplitudes,cur_coil_amplitudes);

            % Check enough repeats
            valid_ops = ~(isnan(cur_observed_amplitudes) | isnan(cur_coil_amplitudes));
            cur_observed_amplitudes = cur_observed_amplitudes(valid_ops);
            cur_coil_amplitudes = cur_coil_amplitudes(valid_ops);
            test_op_idx = op_idx(valid_ops);

            if numel(cur_observed_amplitudes) < 3
                continue
            end

            % Try to remove outliers
            improved = true;
            while current_amplitude_corr < corr_threshold && improved && numel(cur_observed_amplitudes) > 2
                improved = false;
                best_corr = current_amplitude_corr;
                best_remove_idx = [];

                for sel_op_idx = 1:numel(cur_observed_amplitudes)
                    test_observed_amplitudes = cur_observed_amplitudes;
                    test_coil_amplitudes = cur_coil_amplitudes;

                    test_observed_amplitudes(sel_op_idx) = [];
                    test_coil_amplitudes(sel_op_idx) = [];

                    test_corr = corr(test_observed_amplitudes, test_coil_amplitudes);

                    if test_corr > best_corr
                        best_corr = test_corr;
                        best_remove_idx = sel_op_idx;
                    end
                end

                if ~isempty(best_remove_idx)
                    cur_observed_amplitudes(best_remove_idx) = [];
                    cur_coil_amplitudes(best_remove_idx) = [];
                    test_op_idx(best_remove_idx) = [];

                    current_amplitude_corr = best_corr;
                    improved = true;
                end
            end
            grouped_amplitude_corr(sens_idx,axis_idx) = current_amplitude_corr;

            % Remove bad operations
            if current_amplitude_corr < corr_threshold
                phase_relation(sens_idx, axis_idx, :) = NaN;
            else
                removed_ops = setdiff(op_idx, test_op_idx);
                if any(removed_ops)
                    if any(~isnan(phase_relation(sens_idx, axis_idx, removed_ops)))
                        phase_relation(sens_idx, axis_idx, removed_ops) = NaN;
                    end
                end
            end
        end
    end
end
    
% Check coil-wise phase consistency
event_design = S.coil_trigger.event_design;
all_coils = unique(event_design);
for sens_idx = 1:num_sens
    for axis_idx = 1:num_axis
        for coil_idx = 1:length(all_coils)
            op_idx = find(ismember(event_design, all_coils(coil_idx)));
            cur_op_phase_relation = squeeze(phase_relation(sens_idx, axis_idx, op_idx));
            cur_valid_phase_relation = cur_op_phase_relation(~isnan(cur_op_phase_relation));
            if numel(cur_valid_phase_relation) < 2
                continue
            end

            % Check agreement
            if numel(unique(cur_valid_phase_relation)) > 1
                phase_relation(sens_idx, axis_idx, op_idx) = NaN;
            end
        end
    end
end

% Identify similar observations to improve optimiser
use_similarity = S.use_similarity;
if use_similarity
    observed_amplitudes_copy = observed_amplitudes;
    [num_sens, num_axes, ~] = size(observed_amplitudes_copy);
    num_pairs = num_sens * num_axes;
    
    valid_idx = (abs(observed_amplitudes_copy) >= amplitude_range(1)) & ...
                       (abs(observed_amplitudes_copy) <= amplitude_range(2));
    observed_amplitudes_copy(~valid_idx) = NaN;
    
    similarity_matrix = nan(num_pairs, num_pairs);
    for pair_a = 1:num_pairs
        pair_a_sens = mod(pair_a - 1, num_sens) + 1;
        pair_a_axis = ceil(pair_a / num_sens);
    
        for pair_b = pair_a:num_pairs
            pair_b_sens = mod(pair_b - 1, num_sens) + 1;
            pair_b_axis = ceil(pair_b / num_sens);
    
            if pair_a_sens == pair_b_sens
                similarity_matrix(pair_a, pair_b) = NaN;
                continue
            end
    
            phase_pair_a = squeeze(phase_relation(pair_a_sens,pair_a_axis,:));
            phase_pair_b = squeeze(phase_relation(pair_b_sens,pair_b_axis,:));
            amp_pair_a = squeeze(observed_amplitudes_copy(pair_a_sens,pair_a_axis,:));
            amp_pair_b = squeeze(observed_amplitudes_copy(pair_b_sens,pair_b_axis,:));
            valid_idx = ~isnan(phase_pair_a) & ~isnan(phase_pair_b) &...
                        ~isnan(amp_pair_a) & ~isnan(amp_pair_b);
            if sum(valid_idx) < min_valid_for_similar
                similarity_matrix(pair_a, pair_b) = NaN;
                continue
            end
    
            % Phase correlation
            phase_similarity = corr(phase_pair_a(:), phase_pair_b(:), 'rows','complete');
            phase_similarity = abs(phase_similarity);
    
            % Amplitude correlation
            amp_similarity = corr(amp_pair_a(:), amp_pair_b(:), 'rows','complete');
            if isnan(amp_similarity)
                amp_similarity = 0;
            end
    
            % Single measure
            similarity = mean([phase_similarity, amp_similarity]);
            similarity_matrix(pair_a, pair_b) = similarity;
            similarity_matrix(pair_b, pair_a) = similarity;
        end
    end
    
    % Useful order of sensors
    sensor_similarity = nan(num_sens);
    for pair_a_sens = 1:num_sens
        idx_pair_a = pair_a_sens:num_sens:num_pairs;
        for pair_b_sens = 1:num_sens
            if pair_a_sens == pair_b_sens
                continue
            end
            idx_pair_b = (pair_b_sens:num_sens:num_pairs);
            block = similarity_matrix(idx_pair_a, idx_pair_b);
            sensor_similarity(pair_a_sens, pair_b_sens) = sum(block(:), 'omitnan');
        end
    end
    overall_similarity = sum(sensor_similarity, 2, 'omitnan');
    [similarity_vals, sensor_priority] = sort(overall_similarity, 'descend');
    sensor_priority = [sensor_priority(~isnan(similarity_vals)); sensor_priority(isnan(similarity_vals))];
    similarity_available = false(length(similarity_matrix),1);
else
    sensor_priority = (1:num_sens)';
end

% Run estimation
position = cell(num_sens,1);
count = 1;
for sens_idx = sensor_priority'
    fprintf('Processing sensor %d of %d\r', count, num_sens);
    count = count + 1;
    
    start_provided = false;
    if use_similarity
        % Try and find good starting points
        similarity_idx = sens_idx:num_sens:num_pairs;
        num_similarity_idx = length(similarity_idx);
        best_similarity_idx = nan(num_similarity_idx,1);
        for axis_similarity_idx = 1:num_similarity_idx
            similarity_values = similarity_matrix(similarity_idx(axis_similarity_idx),:);
            similarity_values(~similarity_available) = 0;
            similarity_values(similarity_values < good_similarity_threshold) = 0;
            
            % Find maximum valid similarity
            [best_val, best_idx] = max(similarity_values);
            if best_val > 0
                best_similarity_idx(axis_similarity_idx) = best_idx;
            end
        end
    
        % Collect useful positions/orientations
        starting_position = nan(num_axis,3);
        starting_orientation = eye(3,num_axis);
        for axis_similarity_idx = 1:num_similarity_idx
            if ~isnan(best_similarity_idx(axis_similarity_idx))
                cur_similar_sens = mod(best_similarity_idx(axis_similarity_idx) - 1, num_sens) + 1;
                cur_similar_axis = ceil(best_similarity_idx(axis_similarity_idx) / num_sens);
    
                similar_position_table = position{cur_similar_sens};

                similar_sensor_info = sensor_info(cur_similar_sens);
                cur_chan_idx = similar_sensor_info.axis.(axes_labels{cur_similar_axis});
                cur_chan_name = all_chans{cur_chan_idx};

                similar_axis_row_idx = ismember(similar_position_table.name,cur_chan_name);
                similar_position_row = similar_position_table(similar_axis_row_idx,:);
                starting_position(axis_similarity_idx,:) = [similar_position_row.Px, similar_position_row.Py, similar_position_row.Pz];
                starting_orientation(:,axis_similarity_idx) = [similar_position_row.Ox, similar_position_row.Oy, similar_position_row.Oz];
            
                start_provided = true;
            end
        end
        starting_position = mean(starting_position, 1, 'omitnan');
    end
   
    % Select all channels from sensor
    cur_sensor_info = sensor_info(sens_idx);
    cur_chan_indices = nan(1,num_axis);
    for axis_idx = 1:num_axis
        axis_label = axes_labels{axis_idx};
        if isfield(cur_sensor_info.axis, axis_label) && ~isempty(cur_sensor_info.axis.(axis_label))
            cur_chan_indices(axis_idx) = cur_sensor_info.axis.(axis_label);
        end
    end

    % Select only valid axes
    valid_axis = ~isnan(cur_chan_indices);
    enough_valid_ops = sum(~isnan(squeeze(phase_relation(sens_idx,:,:))), 2) >= S.estimation.min_ops;
    valid_axis = find(logical(valid_axis .* enough_valid_ops'));

    if ~isempty(S.highpass)
        sens_dat = proc_D(meg_idx(cur_chan_indices(valid_axis)),:,:);
    else
        sens_dat = S.D(meg_idx(cur_chan_indices(valid_axis)),:,:);
    end
    
    % Plot
    figure(signal_figure_handle)
    for axis_idx = 1:num_axis
        if ~isempty(sub_h{2 + axis_idx})
            cla(sub_h{2 + axis_idx},'reset');
        end
    end
    for valid_axis_idx = 1:length(valid_axis)
        sub_h{valid_axis(valid_axis_idx) + 2} = subplot(2 + num_axis, 1, valid_axis(valid_axis_idx) + 2);

        plot(S.D.time(bounds(1):bounds(2)), sens_dat(valid_axis_idx, bounds(1):bounds(2)),...
        'Color', axis_plotting_colours(rem(valid_axis(valid_axis_idx)-1,size(axis_plotting_colours,1))+1, :))

        y_lims = [min(sens_dat(valid_axis_idx, bounds(1):bounds(2))),max(sens_dat(valid_axis_idx, bounds(1):bounds(2)))];
        ylim(y_lims + [-range(y_lims), range(y_lims)] * 0.1)
        title(axes_labels{valid_axis(valid_axis_idx)})
    end
    drawnow

    % Estimation
    Se = [];
    Se.freq = freq;
    Se.fs = fs;
    Se.amplitude_range = amplitude_range;
    Se.coil_positions = op_coil_positions;
    Se.coil_orientations = op_coil_orientations;
    Se.coil_amplitudes = op_coil_moments;
    Se.phase_relation = squeeze(phase_relation(sens_idx,valid_axis,:));
    Se.observed_amplitudes = squeeze(observed_amplitudes(sens_idx,valid_axis,:));
    Se.plot_output = S.plot_output;
    Se.figure_handle = figure_handle;
    Se.estimation = S.estimation;
    Se.fast = S.fast;
    Se.axes_labs = axes_labels(valid_axis);
    Se.axis_plotting_colours = axis_plotting_colours(valid_axis,:);
    if start_provided && use_similarity
        Se.estimation.starting_pos = starting_position;
        Se.estimation.starting_ori = starting_orientation(:,valid_axis);
    end
    [position{sens_idx}] = estimate_mag_transform(Se);

    % Update similarity and add sensor label back in
    for cur_axis_idx = 1:height(position{sens_idx})
        if use_similarity
            cur_axis_name = position{sens_idx}.name{cur_axis_idx};
            for axis_idx = valid_axis
                if strcmp(axes_labels{axis_idx}, cur_axis_name)
                    axis_match = axis_idx;
                end
            end
            similarity_available(similarity_idx(axis_match)) = true;
        end
        cur_chan_idx = cur_sensor_info.axis.(axes_labels{axis_match});
        cur_chan_name = all_chans{cur_chan_idx};
        position{sens_idx}.name{cur_axis_idx} = cur_chan_name;
    end
end
fprintf('\nAll sensors processed.\n')

positions = vertcat(position{:});

% Set unit from m to mm (in keeping with earlier position file standards)
disp('Setting unit to mm')
positions.Px = positions.Px * 1000;
positions.Py = positions.Py * 1000;
positions.Pz = positions.Pz * 1000;

% Gain is inverse
positions.gain = 1 ./ positions.gain;

if ~isempty(S.output_filename)
    if ~iscell(S.output_filename)
        disp(['Writing positions table to: ',S.output_filename])
        writetable(positions,S.output_filename,'FileType','text','Delimiter','\t');
    else
        disp('Writing positions table to:')
        for file_idx = 1:length(S.output_filename)
            disp(['    ',S.output_filename{file_idx}])
            writetable(positions,S.output_filename{file_idx},'FileType','text','Delimiter','\t');
        end
    end
end

%% Apply calibration
S_c = [];
S_c.D = S.D;
S_c.positions = positions;
S_c.channels = S.channels;
S_c.balance = S.balance;
cal_D = spm_opm_apply_calibration(S_c);
end


%==========================================================================
% - I N P U T   R E A D   O U T
%==========================================================================
function S = input_read_out(S)
% defaults
if ~isfield(S, 'D') || isempty(S.D)
    error('No D object provided.')
end
if ~isfield(S, 'chan_name_format') || isempty(S.chan_name_format)
    S.chan_name_format = '^(X|Y|Z)(.*?)$';
end
if ~isfield(S, 'amplitude_range') || isempty(S.amplitude_range)
    S.amplitude_range = [0.1 1000]*1e3;
end
if ~isfield(S, 'freq') || isempty(S.freq)
    S.freq = 10;
end
if ~isfield(S, 't_delay') || isempty(S.t_delay)
    S.t_delay = 'auto';
end
if (~isfield(S, 'phase_tol') || isempty(S.phase_tol)) && strcmp(S.t_delay, 'auto')
    S.phase_tol = 'auto';
end
if ~isfield(S, 'coil_properties') || isempty(S.coil_properties)
   S.coil_properties.positions = [60, 103.923, 0;
    103.926, 60, 0;
    120, 0, 0;
    103.923, -60, 0;
    60, -103.923, 0;
    0, -120, 0;
    -60, -103.923, 0;
    -103.923, -60, 0;
    -120, 0, 0;
    -103.923, 60, 0;
    -60, 103.923, 0;
    0, 120, 0;
    0, 63.5, 0;
    0, 0, 0;
    -54.99, -31.75, 0;
    54.99, -31.75, 0] ./ 1000;
    S.coil_properties.orientations = repmat([0, 0, 1],[length(S.coil_properties.positions),1]);
    S.coil_properties.resistance = 2120;
    S.coil_properties.inner_r = 3.68/1000;
    S.coil_properties.outer_r = 10.94/1000;
    S.coil_properties.layers = 6;
    S.coil_properties.layer_turns = 45;
end
if ~isfield(S, 'coil_trigger') || isempty(S.coil_trigger)
    S.coil_trigger = struct();
end
if ~isfield(S.coil_trigger, 'signal_channel') || isempty(S.coil_trigger.signal_channel)
    S.coil_trigger.signal_channel = 'A16';
end
if ~isfield(S.coil_trigger, 'event_design') || isempty(S.coil_trigger.event_design)
    S.coil_trigger.event_design = repmat(1:16,[1 4]);
end
if ~isfield(S.coil_trigger, 'event_duration') || isempty(S.coil_trigger.event_duration)
    S.coil_trigger.event_duration = 1;
end
if ~isfield(S.coil_trigger, 'event_trim') || isempty(S.coil_trigger.event_trim)
    S.coil_trigger.event_trim = [1/S.freq, 1/S.freq];
end
if ~isfield(S.coil_trigger, 'event_channel') || isempty(S.coil_trigger.event_channel)
    S.coil_trigger.event_channel = 'A8';
end
if isempty(S.amplitude_range)
    S.amplitude_range = [0 inf];
end
if ~isfield(S, 'plot_output') || isempty(S.plot_output)
    S.plot_output = true;
end
if ~isfield(S, 'balance') || isempty(S.balance)
    S.balance = false;
end
if ~isfield(S, 'highpass')
    S.highpass = 5;
end
if ~isfield(S, 'channels')
    S.channels = [];
end
if ~isfield(S, 'use_similarity') || isempty(S.use_similarity)
    S.use_similarity = true;
end
if ~isfield(S, 'fast') || isempty(S.fast)
    if ~S.use_similarity
        S.fast = false;
    else
        S.fast = true;
    end
end
if ~isfield(S, 'output_filename') || isempty(S.output_filename)
    S.output_filename = [];
end
if ~isfield(S, 'estimation') || isempty(S.estimation)
    S.estimation = struct();
end
default_estimation = struct('permutations', 6,...
                    'min_ops', 8,...
                    'sub_sel_num', 5,...
                    'phase_tol', 'auto',...
                    'pos_bounds', [-0.15, 0.15; -0.15, 0.15; -0.3, -0.05],...
                    'gain_bounds', [0.75, 1.5],...
                    'parallel', false);
fields = fieldnames(default_estimation);
for i = 1:numel(fields)
    fname = fields{i};
    if ~isfield(S.estimation, fname) || isempty(S.estimation.(fname))
        S.estimation.(fname) = default_estimation.(fname);
    end
end


disp('------------------------------------------------------------');
disp('           Calibration Settings Being Used');
disp('------------------------------------------------------------');

disp('General settings:');
disp(['  Channel name format:   [' S.chan_name_format ']']);
disp(['  Amplitude range:       [' num2str(S.amplitude_range(1)) ', ' num2str(S.amplitude_range(2)) ']']);
disp(['  Frequency:             ' num2str(S.freq) ' Hz']);
disp(['  Time delay:            ' num2str(S.t_delay)]);
disp(['  Phase tolerance:       ' num2str(S.phase_tol)]);
disp(['  Coil Trigger channel:  ' S.coil_trigger.event_channel]);
disp(['  Trigger Event Design:  ' num2str(S.coil_trigger.event_design)]);
disp(['  Coil Signal channel:   ' S.coil_trigger.signal_channel]);
disp(['  High-pass frequency:   ' num2str(S.highpass) ' Hz']);
if S.plot_output
    disp('  Plot output:           true');
else
    disp('  Plot output:           false');
end
if S.use_similarity
    disp('  Use similarity:        true');
else
    disp('  Use similarity:        false');
end
if S.fast
    disp('  Fast mode:             true');
else
    disp('  Fast mode:             false');
end
if S.balance
    disp('  Balance:               true');
else
    disp('  Balance:               false');
end
if isempty(S.output_filename)
    disp('  Output filename:       (none)');
else
    if ~iscell(S.output_filename)
        disp(['  Output filename:       ' S.output_filename]);
    else
        disp('  Output files:')
        for file_idx = 1:length(S.output_filename)
            disp(['        ',S.output_filename{file_idx}])
        end
    end
end
disp(' ');
disp('Coil properties:');
disp(['  Number of coils:       ' num2str(size(S.coil_properties.positions,1))]);
disp(['  Resistance:            ' num2str(S.coil_properties.resistance) ' ohm']);
disp(['  Inner radius:          ' num2str(S.coil_properties.inner_r) ' m']);
disp(['  Outer radius:          ' num2str(S.coil_properties.outer_r) ' m']);
disp(['  Layers:                ' num2str(S.coil_properties.layers)]);
disp(['  Turns per layer:       ' num2str(S.coil_properties.layer_turns)]);
disp(' ');
disp('Estimation settings:');
disp(['  Permutations:          ' num2str(S.estimation.permutations)]);
disp(['  Minimum ops:           ' num2str(S.estimation.min_ops)]);
disp(['  Sub-selection count:   ' num2str(S.estimation.sub_sel_num)]);
disp(['  Phase tolerance:       ' num2str(S.estimation.phase_tol)]);
disp(['  Position bounds:       X [' num2str(S.estimation.pos_bounds(1,1)) ', ' ...
                                           num2str(S.estimation.pos_bounds(1,2)) ']  ' ...
                                      'Y [' num2str(S.estimation.pos_bounds(2,1)) ', ' ...
                                           num2str(S.estimation.pos_bounds(2,2)) ']  ' ...
                                      'Z [' num2str(S.estimation.pos_bounds(3,1)) ', ' ...
                                           num2str(S.estimation.pos_bounds(3,2)) ']']);
disp(['  Gain bounds:           [' num2str(S.estimation.gain_bounds(1)) ', ' ...
                                           num2str(S.estimation.gain_bounds(2)) ']']);
if S.estimation.parallel
    disp('  Use parallel:          true');
else
    disp('  Use parallel:          false');
end
disp('------------------------------------------------------------');
disp(' ');
end


%==========================================================================
% - E S T I M A T E   M A G   T R A N S F O R M
%==========================================================================
function [calibration] = estimate_mag_transform(S)
% Estimates the position and orientation of single and multi-axis magnetometer 
% relative to a set of fixed magnetic dipole sources by fitting observed amplitude 
% and phase of the source. Both amplitude and phase errors are used in the objective function.
%
% Inputs:
%   amplitude_range     - 1x2 vector defining valid signal amplitude range [min, max] (V)
%   coil_positions      - n x 3 matrix of dipole (coil) positions (m)
%   coil_amplitudes     - n x 1 vector of relative dipole amplitudes (A·m²)
%   coil_orientations   - n x 3 dipole moment direction unit vector
%   observed_amplitudes - n x m matrix of observed amplitudes (~fT)
%   phase_relation      - n x m matrix of phase/anti-phase conditions
%                           (logical/nan)
%   plot_output       - boolean
%   figure_handle     - figure handle
%   permute           - boolean, permute subselections and check for agreement.
%   fast              - boolean, skip permutations if starting conditions
%                           present
%   use_idx_per_axis  - (optional) internal use
%   iter_lim          - (optional) internal use
%   axes_labs         - n x 1 axes labels
%   estimation: structure with fields:
%       permutations:  number of permutations to run
%       min_ops:       minimum number of valid operations to estimate axis
%       sub_sel_num:   number of operations to estimate for initial estimate/s
%       phase_tol:     tolerance for phase filter
%       pos_bounds:    lower/upper bounds per dimension (m)
%       gain_bounds:   lower/upper bounds for gain (all channels)
%       starting_pos:  1 x 3, position vector
%       starting_ori:  m x 3, orientation vector per axis
%       starting_gain: m x 1, gain per axis
%       parallel:      boolean, whether to use parallel processing where
%                      possible
%
% Outputs:
%       positions : table with fields
%           - name  : sensor axis label
%           - Px, Py, Pz : sensor position (m)
%           - Ox, Oy, Oz : axis orientation (unit vector)
%           - gain       : gain
%           - cal_error  : fitting error
%
% Note:
%   Axes and coil operations with amplitudes outside the specified range
%       are excluded from the fitting.
%   The 'name' field in positions is just the axis - you need to update it
%       to include sensor information at a later stage. 
%   The cal_error is reported per axis, but the optimisation is based on a
%       combined cal_error term for all 'good' axes. 
%   cal_error is the RMS of relative (signed) amplitude error between modelled 
%       and observed responses for the calibration used (position/
%       orientations/gains). If bad channel detection works effectively, no 
%       values higher than ~0.1 should be present, but no threshold is set 
%       in this function.
%
% Author: Nicholas Alexander (n.alexander@ucl.ac.uk)
% Copyright: Department of Imaging Neuroscience, UCL, 2025

% Inputs
min_ops = S.estimation.min_ops;
sub_sel_num = S.estimation.sub_sel_num;
permutations = S.estimation.permutations;
pos_bounds = S.estimation.pos_bounds;
lb_position = pos_bounds(:,1)';
ub_position = pos_bounds(:,2)';
gain_bounds = S.estimation.gain_bounds;
lb_gain = gain_bounds(1);
ub_gain = gain_bounds(2);

amplitude_range = S.amplitude_range;
coil_positions = S.coil_positions;
coil_amplitudes = S.coil_amplitudes;
coil_orientations = S.coil_orientations;
phase_relation = S.phase_relation;
observed_amplitudes = S.observed_amplitudes;
parallel = S.estimation.parallel;
axis_plotting_colours = S.axis_plotting_colours;
axes_labs = S.axes_labs;


% Starting parameters (optional)
num_axes = length(phase_relation(:,1));
default_starting_pos  = mean([lb_position; ub_position]);
default_starting_ori  = eye(3, num_axes);
default_starting_gain = ones(num_axes, 1);

fast = S.fast;
start_provided = false;
if isfield(S.estimation, 'starting_pos') && ~isempty(S.estimation.starting_pos)
    starting_pos = S.estimation.starting_pos;
    start_provided = true;
else
    starting_pos = default_starting_pos;
end
if isfield(S.estimation, 'starting_ori') && ~isempty(S.estimation.starting_ori)
    starting_ori = S.estimation.starting_ori;
    start_provided = true;
else
    starting_ori = default_starting_ori;
end
if isfield(S.estimation, 'starting_gain') && ~isempty(S.estimation.starting_gain)
    starting_gain = S.estimation.starting_gain;
    start_provided = true;
else
    starting_gain = default_starting_gain;
end
if isfield(S, 'iter_lim') && ~isempty(S.iter_lim)
	iter_lim = S.iter_lim;
else
	iter_lim = [0, 5];
end

% Add amplitude filter
if isfield(S, 'use_idx_per_axis') && ~isempty(S.use_idx_per_axis)
    use_idx_per_axis = S.use_idx_per_axis;
else
    use_idx_per_axis = (abs(observed_amplitudes) >= amplitude_range(1)) & ...
                   (abs(observed_amplitudes) <= amplitude_range(2)) & ...
                   ~isnan(phase_relation);
end

axis_use_idx = sum(use_idx_per_axis,2) >= min_ops;
use_idx_per_axis(~axis_use_idx,:) = false;

if all(~axis_use_idx)
    disp('No operations meet the amplitude criteria for any axis.');
    calibration = [];
    return
end

gain_plot_scalar = 0.01;
exp_param_scalar = 1.5;

% Bad axis parameters
improvement_threshold = 3;
relative_difference_threshold = 2;

if ~start_provided || ~fast
    % Global optimisation (permutations)
    ga_options = optimoptions('ga', ...
        'Display', 'none', ...
        'UseParallel', parallel, ...
        'PopulationSize', 256, ...
        'MaxGenerations', 256, ...
        'FunctionTolerance', 1e-2);
    separate_ga_options = ga_options;
    
    bad_axis = ~axis_use_idx;
    first_run = true;
    previous_min_error = inf;
    perm_axis_errors = nan(length(axis_use_idx),permutations);
    while first_run || any(bad_axis)
        if ~first_run
            disp('Bad axis detected. Rerunning without that axis.')
            perm_axis_errors(:,bad_axis) = NaN;
        end
        
        axis_use_idx(bad_axis) = 0;
        num_use_axis = sum(axis_use_idx);
    
        % Setup optimisation
        lb_vector = repmat([-1, -1, -1, lb_gain], [1, num_use_axis]);
        ub_vector = repmat([1,  1,  1, ub_gain], [1, num_use_axis]);
        
        lb = [lb_position, lb_vector];
        ub = [ub_position, ub_vector];
    
        starting_vec = [starting_ori(:,axis_use_idx); starting_gain(axis_use_idx)'];
        starting_param = [starting_pos(:)', starting_vec(:)'];
        n_params = length(starting_param);

        ga_options.InitialPopulationMatrix = starting_param;
        
        % Subselection for permutation
        sub_sel_idx_per_axis_all = cell(permutations,1);
        global_errors = zeros(permutations,1);
        if first_run
            perm_axis_errors = nan(permutations,length(axis_use_idx));
        end
        global_best_all = zeros(permutations,n_params);
        for perm_idx = 1:permutations
            sub_sel_idx_per_axis = false(size(use_idx_per_axis));
            rng(perm_idx);
        
            for axis_idx = 1:length(axis_use_idx)
                if axis_use_idx(axis_idx)
                    valid_ops = find(use_idx_per_axis(axis_idx,:));
                    num_valid = numel(valid_ops);
                    if num_valid <= sub_sel_num
                        sub_sel = valid_ops;
                    else
                        sub_sel = randsample(valid_ops, sub_sel_num);
                    end
                    sub_sel_idx_per_axis(axis_idx, sub_sel) = true;
                end
            end
            sub_sel_idx_per_axis_all{perm_idx} = sub_sel_idx_per_axis;
        
            obj_fun_perm = @(params) objective(params, ...
                coil_positions, coil_amplitudes, coil_orientations, ...
                observed_amplitudes(axis_use_idx,:), ...
                sub_sel_idx_per_axis(axis_use_idx,:), phase_relation(axis_use_idx,:));
        
            % Global for sub selection
            [best_params, best_cost] = ga(obj_fun_perm, n_params, [], [], [], [], lb, ub, [], ga_options);
            global_best_all(perm_idx,:) = best_params;
            global_errors(perm_idx) = best_cost;

            if first_run
                % Per axis error
                for axis_idx = 1:num_axes
                    if ~bad_axis(axis_idx)
                        separate_num_use_axis = 1;
                        separate_axis_idx = false(1,num_axes);
                        separate_axis_idx(axis_idx) = true;
            
                        % Setup optimisation
                        separate_axis_lb_vector = repmat([-1, -1, -1, lb_gain], [1, separate_num_use_axis]);
                        separate_axis_ub_vector = repmat([1,  1,  1, ub_gain], [1, separate_num_use_axis]);
                        
                        separate_axis_lb = [lb_position, separate_axis_lb_vector];
                        separate_axis_ub = [ub_position, separate_axis_ub_vector];
                    
                        separate_axis_starting_vec = [starting_ori(:,separate_axis_idx); starting_gain(separate_axis_idx)'];
                        separate_axis_starting_param = [starting_pos(:)', separate_axis_starting_vec(:)'];
                        separate_axis_n_params = length(separate_axis_starting_param);
                        separate_ga_options.InitialPopulationMatrix = separate_axis_starting_param;
                        
                        % Full selection for local
                        obj_fun_perm = @(params) objective(params, ...
                                    coil_positions, coil_amplitudes, coil_orientations, ...
                                    observed_amplitudes(separate_axis_idx,:), ...
                                    use_idx_per_axis(separate_axis_idx,:), phase_relation(separate_axis_idx,:));
                        [~, perm_axis_errors(perm_idx,axis_idx)] = ga(obj_fun_perm, separate_axis_n_params, [], [], [], [], separate_axis_lb, separate_axis_ub, [], separate_ga_options);
                    end
                end
            end
        end
        first_run = false;

%        % Variability across permutations - worth implementing as bad channel
%        % check?
%         norm_global_best_all = (global_best_all - mean(global_best_all,1) ) ./ mean(global_best_all,1);

        % Check improvement
        [current_min_error, best_perm_idx] = min(global_errors);
        global_best = global_best_all(best_perm_idx,:);
        if ~isinf(previous_min_error)
            improvement = previous_min_error / current_min_error;
            if improvement < improvement_threshold
                disp('No improvement after removing that axis. Reverting back.');
                axis_use_idx = previous_axis_use_idx;
                global_best = previous_optimal_params;
                break
            else
                if sum(~bad_axis) == 1
                    break
                end
            end
        end
        previous_min_error = current_min_error;
        previous_axis_use_idx = axis_use_idx;
        previous_optimal_params = global_best;

        % Check for bad axes
        best_cost_axis = min(perm_axis_errors);
        [max_axis_val, max_axis_idx] = max(best_cost_axis);
        others = best_cost_axis(best_cost_axis ~= max_axis_val);
        if all(max_axis_val > relative_difference_threshold * others)
            bad_axis(max_axis_idx) = true;
        end
    end

    num_use_axis = sum(axis_use_idx);
    lb_vector = repmat([-1, -1, -1, lb_gain], [1, num_use_axis]);
    ub_vector = repmat([1,  1,  1, ub_gain], [1, num_use_axis]);
    lb = [lb_position, lb_vector];
    ub = [ub_position, ub_vector];
    

    % Local
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'none', 'MaxFunctionEvaluations', 1e8, ...
                'MaxIterations', 1e8, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
    
    % Full selection for local
    obj_fun = @(params) objective(params, ...
                coil_positions, coil_amplitudes, coil_orientations, ...
                observed_amplitudes(axis_use_idx,:), ...
                use_idx_per_axis(axis_use_idx,:), phase_relation(axis_use_idx,:));
    [optimal_params] = fmincon(obj_fun, global_best, [], [], [], [], lb, ub, [], options);
else
    bad_axis = ~axis_use_idx;
    first_run = true;
    previous_min_error = inf;
    previous_optimal_params = [];
    while first_run || any(bad_axis)
        if ~first_run
            disp('Bad axis detected. Rerunning without that axis.')
        end
        first_run = false;
    
        axis_use_idx(bad_axis) = 0;
        num_use_axis = sum(axis_use_idx);
        
        % Setup optimisation
        lb_vector = repmat([-1, -1, -1, lb_gain], [1, num_use_axis]);
        ub_vector = repmat([1,  1,  1, ub_gain], [1, num_use_axis]);
        
        lb = [lb_position, lb_vector];
        ub = [ub_position, ub_vector];
        
        starting_vec = [starting_ori(:,axis_use_idx); starting_gain(axis_use_idx)'];
        starting_param = [starting_pos(:)', starting_vec(:)'];

        options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'none', 'MaxFunctionEvaluations', 1e8, ...
                    'MaxIterations', 1e8, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
        
        % Full selection for local
        obj_fun = @(params) objective(params, ...
                    coil_positions, coil_amplitudes, coil_orientations, ...
                    observed_amplitudes(axis_use_idx,:), ...
                    use_idx_per_axis(axis_use_idx,:), phase_relation(axis_use_idx,:));
        [optimal_params, best_cost] = fmincon(obj_fun, starting_param, [], [], [], [], lb, ub, [], options);

        best_cost_axis = nan(num_axes,1);
        for axis_idx = 1:num_axes
            if ~bad_axis(axis_idx)
                separate_num_use_axis = 1;
                separate_axis_idx = false(1,num_axes);
                separate_axis_idx(axis_idx) = true;
    
                % Setup optimisation
                lb_vector = repmat([-1, -1, -1, lb_gain], [1, separate_num_use_axis]);
                ub_vector = repmat([1,  1,  1, ub_gain], [1, separate_num_use_axis]);
                
                lb = [lb_position, lb_vector];
                ub = [ub_position, ub_vector];
            
                starting_vec = [starting_ori(:,separate_axis_idx); starting_gain(separate_axis_idx)'];
                starting_param = [starting_pos(:)', starting_vec(:)'];
    
                
                options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'none', 'MaxFunctionEvaluations', 1e8, ...
                            'MaxIterations', 1e8, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
                
                % Full selection for local
                obj_fun = @(params) objective(params, ...
                            coil_positions, coil_amplitudes, coil_orientations, ...
                            observed_amplitudes(separate_axis_idx,:), ...
                            use_idx_per_axis(separate_axis_idx,:), phase_relation(separate_axis_idx,:));
                [~, best_cost_axis(axis_idx)] = fmincon(obj_fun, starting_param, [], [], [], [], lb, ub, [], options);
            end
        end

        % Check improvement
        current_min_error = best_cost;
        if ~isinf(previous_min_error)
            improvement = previous_min_error / current_min_error;
            if improvement < improvement_threshold
                disp('No improvement after removing that axis. Reverting back.');
                axis_use_idx = previous_axis_use_idx;
                optimal_params = previous_optimal_params;
                break
            else
                if sum(~bad_axis) == 1
                    break
                end
            end
        end
        previous_min_error = current_min_error;
        previous_axis_use_idx = axis_use_idx;
        previous_optimal_params = optimal_params;
        
        % Check for bad axes
        [max_axis_val, max_axis_idx] = max(best_cost_axis);
        others = best_cost_axis(best_cost_axis ~= max_axis_val);
        if all(max_axis_val > relative_difference_threshold * others)
            bad_axis(max_axis_idx) = true;
        end
    end
end

% Per axis error calc
[~, axis_errors] = objective(optimal_params, ...
            coil_positions, coil_amplitudes, coil_orientations, ...
            observed_amplitudes(axis_use_idx,:), ...
            use_idx_per_axis(axis_use_idx,:), phase_relation(axis_use_idx,:));

% Unpack
position = optimal_params(1:3);
param_vectors = optimal_params(4:end);
num_axes = length(param_vectors)/4;
orientations = zeros(num_axes,3);
gains = zeros(num_axes,1);
for axis_idx = 1:num_axes
    ori(1) = param_vectors((axis_idx-1)*4 + 1);
    ori(2) = param_vectors((axis_idx-1)*4 + 2);
    ori(3) = param_vectors((axis_idx-1)*4 + 3);
    gain  = param_vectors((axis_idx-1)*4 + 4);

    orientations(axis_idx,:) = ori ./ norm(ori);
    gains(axis_idx) = gain;
end

% Format positions output
cur_axes_labs = axes_labs(axis_use_idx);
n_axes = sum(axis_use_idx);
rows = cell(n_axes, 1);
for axis_idx = 1:n_axes
    name = cur_axes_labs{axis_idx};
    Px = position(1);
    Py = position(2);
    Pz = position(3);
    Ox = orientations(axis_idx,1);
    Oy = orientations(axis_idx,2);
    Oz = orientations(axis_idx,3);
    gain = gains(axis_idx);
    cal_error = axis_errors(axis_idx);
    rows{axis_idx} = table({name}, Px, Py, Pz, Ox, Oy, Oz, gain, cal_error, ...
        'VariableNames', {'name','Px','Py','Pz','Ox','Oy','Oz','gain','cal_error'});
end

% Concatenate all rows into one table
calibration = vertcat(rows{:});

% Try improving estimation
if iter_lim(1) < iter_lim(2)
	iter_lim(1) = iter_lim(1) + 1;

	% Check for amplitude threshold changes due to change in gain
	new_use_idx_per_axis = use_idx_per_axis;
	new_gains = starting_gain;
	count = 0;
	for axis_idx = 1:length(new_gains)
    	if axis_use_idx(axis_idx)
        	count = count + 1;
        	new_gains(axis_idx) = gains(count);
    	end
    	new_amplitude_range = amplitude_range * new_gains(axis_idx);
    	new_use_idx_per_axis(axis_idx,:) = (abs(observed_amplitudes(axis_idx,:)) >= new_amplitude_range(1)) & ...
               	(abs(observed_amplitudes(axis_idx,:)) <= new_amplitude_range(2)) & ...
               	~isnan(phase_relation(axis_idx,:));
	end
	update_amplitude_threshold = ~isequal(new_use_idx_per_axis, use_idx_per_axis);
	if update_amplitude_threshold
    	disp('Gain change of amplitude theshold. Increasing search space.')
    	new_gain_bounds = gain_bounds .* [1/exp_param_scalar, exp_param_scalar];
	else
    	new_gain_bounds = gain_bounds;
	end
	
	% Also check for positions on the boundary
	pos_on_boundary = any([ismember(Px, [lb_position(1), ub_position(1)]), ...
    	ismember(Py, [lb_position(2), ub_position(2)]), ...
    	ismember(Pz, [lb_position(3), ub_position(3)])]);
	if pos_on_boundary
    	disp('Sensor on boundary. Increasing search space.')
    	new_pos_bounds = pos_bounds .* exp_param_scalar;
	else
    	new_pos_bounds = pos_bounds;
	end
	
	% Rerun with updated estimation parameters if necessary
	if pos_on_boundary || update_amplitude_threshold
    	disp('Running new estimation')
    	S_new = S;
    	S_new.use_idx_per_axis = new_use_idx_per_axis;
    	S_new.estimation.starting_gain = new_gains;
    	S_new.estimation.gain_bounds = new_gain_bounds;
    	S_new.estimation.pos_bounds = new_pos_bounds;
    	S_new.plot_output = false;
		S_new.iter_lim = iter_lim;
    	[calibration] = estimate_mag_transform(S_new);
	end
else
	disp('Improvement iteration limit exceeded.')
end


% Plot
if S.plot_output
    if ~isfield(S,'figure_handle') || ~isvalid(S.figure_handle)
        S.figure_handle = figure('Name','HALO Sensor Plot');
    else
        figure(S.figure_handle);
    end
    hold on

    % Plot sensor
    if ~isempty(axis_use_idx) && any(axis_use_idx) && ~isempty(orientations)
        for axis_idx = 1:size(orientations, 1)
            quiver3(position(1), position(2), position(3), ...
                    orientations(axis_idx, 1), ...
                    orientations(axis_idx, 2), ...
                    orientations(axis_idx, 3), ...
                    gain_plot_scalar*gains(axis_idx), ...
                    'color', axis_plotting_colours(axis_idx,:));
        end
        axis equal

        xlim([lb_position(1) ub_position(1)])
        ylim([lb_position(2) ub_position(2)])
        zlim([lb_position(3) ub_position(3)])
        drawnow
    end
end
end

%==========================================================================
% - O B J E C T I V E
%==========================================================================
function [error, axis_errors] = objective(params, ...
    coil_positions, coil_amplitudes, coil_orientations, ...
    observed_amplitudes, use_idx_per_axis, phase_relation)

    mag_position = params(1:3);
    param_vectors = params(4:end);
    num_axes = length(param_vectors)/4;
    mag_orientations = zeros(num_axes,3);
    mag_gains = zeros(num_axes,1);
    
    for axis_idx = 1:num_axes
        ori(1) = param_vectors((axis_idx-1)*4 + 1);
        ori(2) = param_vectors((axis_idx-1)*4 + 2);
        ori(3) = param_vectors((axis_idx-1)*4 + 3);
        gain  = param_vectors((axis_idx-1)*4 + 4);
    
        if norm(ori) < 1e-12
            ori = [1 0 0];
        end
        mag_orientations(axis_idx,:) = ori ./ norm(ori);
        mag_gains(axis_idx) = gain;
    end

    % fT
    mag_scale = 1e15; 

    % Simulate for all operations, per axis
    num_ops = length(coil_positions(:,1));
    V_diff = nan(num_axes,num_ops);
    for axis_idx = 1:num_axes
        % Apply amp/phase filter per axis
        use_idx = use_idx_per_axis(axis_idx, :);
        if ~any(use_idx)
            continue
        end

        % Observation info
        axis_observed_amplitudes = observed_amplitudes(axis_idx, use_idx)';
        axis_mag_gain = mag_gains(axis_idx);
        axis_mag_orientation = mag_orientations(axis_idx,:);
        axis_mag_orientation_scaled = axis_mag_orientation .* axis_mag_gain;

        V = model_sensor_measurement_of_dipole(mag_position, axis_mag_orientation_scaled .* mag_scale, coil_positions(use_idx,:), ...
                    coil_orientations(use_idx,:) .* coil_amplitudes(use_idx));

        axis_phase_relation = phase_relation(axis_idx, use_idx)';
        amp_error = abs(V - (axis_observed_amplitudes .* axis_phase_relation)) ./ axis_observed_amplitudes;

        V_diff(axis_idx,use_idx) = amp_error.^2;
    end
 
%     % Group by coil position - this seems sensible, but needs validation
%     % Use variability within coil position? Also, take most of this
%     % outside the objective for speed. 
%     [~,~,position_group] = unique(coil_positions,'rows');
%     num_unique_pos = max(position_group);
%     V_diff_col = nan(num_axes, num_unique_pos);
% 
%     for pos_idx = 1:num_unique_pos
%         idx = (position_group == pos_idx);
%         for axis_idx = 1:num_axes
%             vals = V_diff(axis_idx, idx);
%             V_diff_col(axis_idx, pos_idx) = mean(vals(~isnan(vals)));
%         end
%     end
%     error = sqrt(mean(V_diff_col(~isnan(V_diff_col))));

    error = sqrt(mean(V_diff(~isnan(V_diff))));

    if nargout > 1
        axis_errors = nan(num_axes, 1);
        for axis_idx = 1:num_axes
            vals = V_diff(axis_idx, :);
            axis_errors(axis_idx) = sqrt(mean(vals(~isnan(vals))));
        end
    end
end

%==========================================================================
% - M O D E L   S E N S O R   M E A S U R E M E N T   O F   D I P O L E
%==========================================================================
function V = model_sensor_measurement_of_dipole(r, o, p, m)
%   r = 1 x 3 sensor position (m)
%   o = 1 x 3 sensor axis unit vector * gain
%   p = n x 3 dipole position (m)
%   m = n x 3 dipole moment (A·m^2)

% Distance from coil origin to sensor position
r_rel = r - p;
delta = sqrt(sum(r_rel.^2, 2));

% Separate gain and orientation
g = norm(o);
ori = o ./ g;
    
% Magnetic dipole field equation
mu0 = 4 * pi * 1e-7;
dot_m_rp = sum(m .* r_rel, 2);

term1 = 3 * r_rel .* dot_m_rp ./ delta.^5;
term2 = m ./ delta.^3;
term3 = mu0 / (4 * pi);
V = g * term3 * sum((term1 - term2) .* ori, 2);
end
