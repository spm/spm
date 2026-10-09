function cfg = dts_meeg_test_defaults(varargin)
% Defaults and path conventions for the SPM DCM M/EEG regression suite.
% Authored by Pranay Yadav in 2026

cfg = struct();
cfg.models = {'ERP','SEP','LFP','CMC','CMM','NMM','MFM','NMDA','CMM_NMDA'};
cfg.integrators = {'spm_int_L','spm_int_B','spm_int_J'};

% NLSI iterations for fixed data cases that need strict estimation.
cfg.strict_fit_Nmax = 32;

% NLSI iterations for broad regression fits that only check the inversion path.
cfg.basic_fit_Nmax = 4;

% FixedDataCase names run in strict mode for all models.
cfg.strict_fixed_data_cases = {'simulated_evoked_lfp', 'sample_lfp'};

cfg.integrator = 'spm_int_L';
cfg.IMG_Nm = 6;
cfg.IMG_radius = 16;
cfg.ECD_location = 0;
cfg.keep_outputs = false;
cfg.preserve_on_failure = true;
cfg.parallel = false;
cfg.pool_workers = 8;
cfg.deterministic_work_root = false;

cfg.thresholds = [];

% Stability bound: largest real Jacobian eigenvalue must stay below this.
cfg.thresholds.jacobian_max_real_eig = 1e-9;

% Maximum normalised deviation from the initialized baseline trajectory.
cfg.thresholds.baseline_error = 0.1;

% Minimum correlation with the spm_int_J forward response.
cfg.thresholds.integrator_min_corr = 0.75;

% Optional R2 floor for self-gen LFP fits; 0 disables it.
cfg.thresholds.self_lfp_min_r2 = 0;

% Allowed absolute free-energy drift in log-evidence units.
cfg.thresholds.F_margin = 5;

% Allowed maximum absolute posterior-parameter drift.
cfg.thresholds.Ep_margin = 0.1;

% Allowed absolute R2 percentage-point drift from the reference fit.
cfg.thresholds.R2_margin = 5;

% Optional absolute R2 floor for strict fits; 0 disables it.
cfg.thresholds.strict_min_R2 = 0;

% Optional absolute R2 floor for basic fits; 0 disables it.
cfg.thresholds.basic_min_R2 = 0;

cfg.suite_root = fullfile(spm('Dir'), 'tests', 'data', 'MEEG', 'DCM_Test_Suite');
cfg.prep_root = fullfile(cfg.suite_root, 'prep');
cfg.ref_root = fullfile(cfg.suite_root, 'refdata');
cfg.reference_file = fullfile(cfg.ref_root, 'dts_meeg_reference.mat');
cfg.output_root = fullfile(cfg.suite_root, 'output', 'tmpdir');

if nargin > 0
    cfg = apply_options(cfg, varargin{:});
end
end

function cfg = apply_options(cfg, varargin)
for i = 1:numel(varargin)
    options = varargin{i};
    if isempty(options), continue, end
    names = fieldnames(options);
    for n = 1:numel(names)
        name = names{n};
        if isfield(cfg, name) && isstruct(cfg.(name)) && isstruct(options.(name)) && ...
                isscalar(cfg.(name)) && isscalar(options.(name))
            nested = fieldnames(options.(name));
            for j = 1:numel(nested)
                cfg.(name).(nested{j}) = options.(name).(nested{j});
            end
        else
            cfg.(name) = options.(name);
        end
    end
end
end
