function cases = dts_meeg_fixed_data_cases(name)
% Fixed prepared datasets used by DCM M/EEG fit regression tests.
% Authored by Pranay Yadav in 2026

cases = [];

cases(1).name = 'simulated_evoked_lfp';
cases(1).test_id = 'test_05_simulated_evoked_lfp';
cases(1).prep_dir = 'test_05_simulated_evoked_lfp';
cases(1).file_name = 'simulated_evoked_lfp.mat';
cases(1).fixture_label = 'simulated_ERP';
cases(1).spatial = 'LFP';
cases(1).generator_model = 'ERP';

cases(2).name = 'synthetic_triphasic_lfp';
cases(2).test_id = 'test_05_synthetic_triphasic_lfp';
cases(2).prep_dir = 'test_05_synthetic_triphasic_lfp';
cases(2).file_name = 'synthetic_triphasic_lfp.mat';
cases(2).fixture_label = 'triphasic';
cases(2).spatial = 'LFP';
cases(2).generator_model = 'triphasic';

cases(3).name = 'sample_lfp';
cases(3).test_id = 'test_05_sample_lfp';
cases(3).prep_dir = 'test_05_sample_lfp';
cases(3).file_name = 'A1_LFP.mat';
cases(3).fixture_label = 'sample_A1_LFP';
cases(3).spatial = 'LFP';
cases(3).generator_model = 'sample_A1_LFP';

cases(4).name = 'sample_sensor_img';
cases(4).test_id = 'test_05_sample_sensor_img';
cases(4).prep_dir = 'local/test_05_sample_sensor_img';
cases(4).file_name = 'EEG_MMN.mat';
cases(4).fixture_label = 'sample_EEG_MMN_IMG';
cases(4).spatial = 'IMG';
cases(4).generator_model = 'sample_EEG_MMN';

cases(5).name = 'sample_sensor_ecd';
cases(5).test_id = 'test_05_sample_sensor_ecd';
cases(5).prep_dir = 'local/test_05_sample_sensor_ecd';
cases(5).file_name = 'EEG_MMN.mat';
cases(5).fixture_label = 'sample_EEG_MMN_ECD';
cases(5).spatial = 'ECD';
cases(5).generator_model = 'sample_EEG_MMN';

if nargin < 1 || isempty(name)
    return
end

idx = find(strcmp({cases.name}, name), 1);
if isempty(idx)
    error('Unknown fixed DCM M/EEG data case: %s', name);
end
cases = cases(idx);
end
