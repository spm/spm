function manifest = dts_prep01_make_test_fixtures(options)
% Generate immutable DCM M/EEG prep fixtures under SPM tests/data.
% Authored by Pranay Yadav in 2026
%
% This is an explicit maintenance script. Runtime tests do not call it.

if nargin < 1, options = struct(); end
support_dir = fileparts(mfilename('fullpath'));
addpath(support_dir);
dts_meeg_setup_paths();
cfg = dts_meeg_test_defaults(options);
cases = dts_meeg_fixed_data_cases();

if ~exist(cfg.prep_root, 'dir'), mkdir(cfg.prep_root); end

manifest = struct();
manifest.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
manifest.spm_dir = spm('Dir');
manifest.prep_root = cfg.prep_root;
manifest.outputs = struct();
manifest.fixed_data_cases = cases;

sim_case = dts_meeg_fixed_data_cases('simulated_evoked_lfp');
sim_dir = fullfile(cfg.prep_root, sim_case.prep_dir);
if ~exist(sim_dir, 'dir'), mkdir(sim_dir); end
y = dts_meeg_prior_lfp_response('ERP', 0.001, 500, cfg.integrator);
manifest.outputs.(sim_case.name) = write_noisy_lfp( ...
    y, fullfile(sim_dir, sim_case.file_name), sim_case.fixture_label);

trip_case = dts_meeg_fixed_data_cases('synthetic_triphasic_lfp');
trip_dir = fullfile(cfg.prep_root, trip_case.prep_dir);
if ~exist(trip_dir, 'dir'), mkdir(trip_dir); end
t = (0:0.001:(0.5 - 0.001))';
g = @(mu, sd, amp) amp*exp(-0.5*((t - mu)/sd).^2);
y = g(0.070, 0.018, 1.0) - g(0.125, 0.025, 1.4) + g(0.205, 0.045, 0.7);
manifest.outputs.(trip_case.name) = write_noisy_lfp( ...
    y, fullfile(trip_dir, trip_case.file_name), trip_case.fixture_label);

sensor_dir = fullfile(cfg.prep_root, 'test_05_sample_sensor_data');
if exist(sensor_dir, 'dir'), rmdir(sensor_dir, 's'); end
mkdir(sensor_dir);

raw_dir = fullfile(spm('Dir'), 'tests', 'data', 'eeg_mmn');
raw_bdf = fullfile(raw_dir, 'subject1.bdf');
raw_pol = fullfile(raw_dir, 'sensors.pol');
if exist(raw_bdf, 'file') ~= 2 || exist(raw_pol, 'file') ~= 2
    error('Missing raw MMN test data. Expected:\n%s\n%s', raw_bdf, raw_pol);
end
eeg_work_dir = fullfile(cfg.output_root, 'prep_work_eeg');
if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
eeg_outputs = dts_meeg_prepare_mmn_eeg(raw_bdf, raw_pol, sensor_dir, eeg_work_dir, true, true, true);

lfp_case = dts_meeg_fixed_data_cases('sample_lfp');
lfp_dir = fullfile(cfg.prep_root, lfp_case.prep_dir);
lfp_file = fullfile(lfp_dir, lfp_case.file_name);
if ~exist(lfp_dir, 'dir'), mkdir(lfp_dir); end
lfp_work_dir = fullfile(cfg.output_root, 'prep_work_lfp');
manifest.outputs.(lfp_case.name) = dts_meeg_prepare_a1_lfp( ...
    eeg_outputs.final_mat, lfp_file, lfp_work_dir, cfg);

D = spm_eeg_load(eeg_outputs.final_mat);
try
    gainmat_file = fullfile(D.path, D.inv{D.val}.gainmat);
catch
    gainmat_file = '';
end
try
    D.inv = {};
catch
end
save(D);

if ~isempty(gainmat_file) && exist(gainmat_file, 'file') == 2
    delete(gainmat_file);
end
extra_gainmats = dir(fullfile(sensor_dir, 'SPMgainmatrix_*.mat'));
for i = 1:numel(extra_gainmats)
    delete(fullfile(extra_gainmats(i).folder, extra_gainmats(i).name));
end
eeg_outputs.gainmat_file = '';
manifest.outputs.sample_sensor_data = eeg_outputs;

save(fullfile(cfg.prep_root, 'prep_manifest.mat'), 'manifest', spm_get_defaults('mat.format'));
local_prep_root = fullfile(cfg.prep_root, 'local');
if exist(local_prep_root, 'dir')
    rmdir(local_prep_root, 's');
end
if exist(cfg.output_root, 'dir')
    rmdir(cfg.output_root, 's');
end
fprintf('\nPrepared DCM M/EEG regression fixtures under:\n%s\n', cfg.prep_root);
end

function data_file = write_noisy_lfp(y, data_file, label)
y = y(:) - y(1);
y = y/std(y);
state = rng;
cleanup = onCleanup(@()rng(state));
rng(1, 'twister');
y = y + randn(size(y))*exp(-6/2);
D = dts_meeg_write_lfp_fixture(y, 0.001, data_file, label);
data_file = D.fullfile();
try
    spm_eeg_load(data_file);
catch ME
    error('Generated M/EEG fixture is not loadable: %s\n%s', data_file, ME.message);
end
end
