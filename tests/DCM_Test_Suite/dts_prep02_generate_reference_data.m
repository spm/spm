function reference = dts_prep02_generate_reference_data(parallel, cleanup_generated)
% Generate deterministic DCM fit references under SPM tests/data.
% Authored by Pranay Yadav in 2026
%
% This is an explicit maintenance script. Runtime tests do not call it.
%
% Input:
%   parallel - optional logical scalar. Default false. When true,
%              model fits use the suite default 8-worker pool.
%   cleanup_generated - optional logical scalar. Default false. When true,
%                       removes generated fit outputs and prep/local after
%                       writing dts_meeg_reference.mat/.csv.

if nargin < 1 || isempty(parallel), parallel = false; end
if nargin < 2 || isempty(cleanup_generated), cleanup_generated = false; end
parallel = logical(parallel);
cleanup_generated = logical(cleanup_generated);
support_dir = fileparts(mfilename('fullpath'));
addpath(support_dir);
dts_meeg_setup_paths();
cfg = dts_meeg_test_defaults(test_regress_spm_dcm_meeg.Options);
cfg.parallel = parallel;

if parallel
    if exist('gcp', 'file') ~= 2 || exist('parpool', 'file') ~= 2
        error('Parallel reference generation requested, but Parallel Computing Toolbox functions are unavailable.');
    end
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(cfg.pool_workers);
    elseif pool.NumWorkers ~= cfg.pool_workers
        delete(pool);
        parpool(cfg.pool_workers);
    end
end

if ~exist(cfg.ref_root, 'dir'), mkdir(cfg.ref_root); end
if ~exist(cfg.prep_root, 'dir')
    error('Missing prep root: %s. Run tests/DCM_Test_Suite/dts_prep01_make_test_fixtures.m first.', cfg.prep_root);
end
dts_meeg_sample_sensor_fixture(cfg, 'genrefdata');

output_root = fullfile(cfg.ref_root, 'fits');
if exist(output_root, 'dir')
    rmdir(output_root, 's');
end
mkdir(output_root);

cases = dts_meeg_fixed_data_cases();
jobs = repmat(struct('data_case', cases(1), 'model', '', 'data_file', ''), ...
    1, numel(cases)*numel(cfg.models));
job_idx = 0;
for c = 1:numel(cases)
    data_file = fullfile(cfg.prep_root, cases(c).prep_dir, cases(c).file_name);
    if exist(data_file, 'file') ~= 2
        error('Missing prepared fixture for %s: %s', cases(c).name, data_file);
    end
    for m = 1:numel(cfg.models)
        job_idx = job_idx + 1;
        jobs(job_idx) = struct( ...
            'data_case', cases(c), ...
            'model', cfg.models{m}, ...
            'data_file', data_file);
    end
end

result_template = struct('test_id', '', 'model', '', 'spatial', '', ...
    'generator_model', '', 'F', NaN, 'R2', NaN, 'Ep', [], ...
    'Nmax', NaN, 'estimation_mode', '', 'min_R2', NaN, ...
    'file', '', 'output_dir', '', 'ok', false, 'error', '');
fixed_results = repmat(result_template, 1, numel(jobs));

if parallel
    parfor i = 1:numel(jobs)
        data_case = jobs(i).data_case;
        fit_cfg = estimation_cfg_for_case(cfg, data_case.name);
        fixed_results(i) = dts_meeg_fit_erp_model(data_case.test_id, jobs(i).model, ...
            jobs(i).data_file, data_case.spatial, data_case.generator_model, output_root, fit_cfg);
    end
else
    for i = 1:numel(jobs)
        data_case = jobs(i).data_case;
        fit_cfg = estimation_cfg_for_case(cfg, data_case.name);
        fixed_results(i) = dts_meeg_fit_erp_model(data_case.test_id, jobs(i).model, ...
            jobs(i).data_file, data_case.spatial, data_case.generator_model, output_root, fit_cfg);
    end
end

models = cfg.models;
csd_results = repmat(result_template, 1, numel(models));
if parallel
    parfor i = 1:numel(models)
        model = models{i};
        fit_cfg = estimation_cfg_for_case(cfg, '');
        csd_results(i) = dts_meeg_fit_csd_model(model, output_root, fit_cfg);
    end
else
    for i = 1:numel(models)
        fit_cfg = estimation_cfg_for_case(cfg, '');
        csd_results(i) = dts_meeg_fit_csd_model(models{i}, output_root, fit_cfg);
    end
end

all_results = [fixed_results csd_results];

for i = 1:numel(all_results)
    r = all_results(i);
    if ~isempty(r.error)
        error('Reference generation failed for %s/%s:\n%s', r.test_id, r.model, r.error);
    elseif ~r.ok
        error('%s/%s did not produce finite F/R2/Ep.', r.test_id, r.model);
    end
end

reference = struct();
reference.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
reference.spm_dir = spm('Dir');
reference.matlab_version = version;
reference.platform = computer;
reference.models = cfg.models;
reference.fixed_data_cases = {cases.name};
reference.csd_data = 'LFP_BG';
reference.test_options = test_regress_spm_dcm_meeg.Options;
reference.thresholds = cfg.thresholds;
reference.rows = repmat(struct( ...
    'test_id', '', 'model', '', 'spatial', '', 'generator_model', '', ...
    'F', NaN, 'R2', NaN, 'Ep', [], 'Nmax', NaN, 'estimation_mode', '', 'file', '', ...
    'F_margin', NaN, 'Ep_margin', NaN, 'R2_abs_tol', NaN, ...
    'min_R2', NaN), 1, numel(all_results));

for i = 1:numel(all_results)
    r = all_results(i);
    row = reference.rows(i);
    row.test_id = r.test_id;
    row.model = r.model;
    row.spatial = r.spatial;
    row.generator_model = r.generator_model;
    row.F = double(full(r.F));
    row.R2 = double(full(r.R2));
    row.Ep = r.Ep;
    row.Nmax = r.Nmax;
    row.estimation_mode = r.estimation_mode;
    row.file = r.file;
    row.F_margin = cfg.thresholds.F_margin;
    row.Ep_margin = cfg.thresholds.Ep_margin;
    row.R2_abs_tol = cfg.thresholds.R2_margin;
    row.min_R2 = r.min_R2;
    reference.rows(i) = row;
end

save(cfg.reference_file, 'reference', spm_get_defaults('mat.format'));
csv_file = fullfile(cfg.ref_root, 'dts_meeg_reference.csv');
fid = fopen(csv_file, 'w');
cleanup = onCleanup(@()fclose(fid));
fprintf(fid, 'test_id,model,spatial,generator_model,estimation_mode,Nmax,F,R2,F_margin,Ep_margin,R2_abs_tol,min_R2,DCM_file\n');
for i = 1:numel(reference.rows)
    r = reference.rows(i);
    fprintf(fid, '%s,%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n', ...
        r.test_id, r.model, r.spatial, r.generator_model, r.estimation_mode, double(full(r.Nmax)), ...
        double(full(r.F)), double(full(r.R2)), double(full(r.F_margin)), ...
        double(full(r.Ep_margin)), double(full(r.R2_abs_tol)), ...
        double(full(r.min_R2)), r.file);
end
clear cleanup

fprintf('\nGenerated DCM M/EEG references under:\n%s\n', cfg.ref_root);

if cleanup_generated
    if exist(output_root, 'dir')
        rmdir(output_root, 's');
    end
    local_prep_root = fullfile(cfg.prep_root, 'local');
    if exist(local_prep_root, 'dir')
        rmdir(local_prep_root, 's');
    end
    fprintf('Removed generated reference fit outputs and prep/local.\n');
end
end

function fit_cfg = estimation_cfg_for_case(cfg, fixed_data_case)
fit_cfg = cfg;
if any(strcmp(fit_cfg.strict_fixed_data_cases, fixed_data_case))
    fit_cfg.estimation_mode = 'strict';
    fit_cfg.Nmax = fit_cfg.strict_fit_Nmax;
    fit_cfg.thresholds.min_R2 = fit_cfg.thresholds.strict_min_R2;
else
    fit_cfg.estimation_mode = 'basic';
    fit_cfg.Nmax = fit_cfg.basic_fit_Nmax;
    fit_cfg.thresholds.min_R2 = fit_cfg.thresholds.basic_min_R2;
end
end
