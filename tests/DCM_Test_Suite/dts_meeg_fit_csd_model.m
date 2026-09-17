function result = dts_meeg_fit_csd_model(model, work_root, cfg)
% Fit one CSD DCM model to SPM's shipped LFP_BG dataset.
% Authored by Pranay Yadav in 2026

test_id = 'test_07_csd_lfp_fit';
estimation_mode = 'unspecified';
if isfield(cfg, 'estimation_mode')
    estimation_mode = cfg.estimation_mode;
end

result = struct( ...
    'test_id', test_id, ...
    'model', model, ...
    'spatial', 'LFP', ...
    'generator_model', 'LFP_BG', ...
    'F', NaN, ...
    'R2', NaN, ...
    'Ep', [], ...
    'Nmax', cfg.Nmax, ...
    'estimation_mode', estimation_mode, ...
    'min_R2', cfg.thresholds.min_R2, ...
    'file', '', ...
    'output_dir', '', ...
    'ok', false, ...
    'error', '');

try
    previous_num_threads = maxNumCompThreads();
    cleanup_threads = onCleanup(@()maxNumCompThreads(previous_num_threads));
    maxNumCompThreads(1);

    data_file = fullfile(spm('Dir'), 'tests', 'data', 'MEEG', 'LFP_BG.mat');
    if exist(data_file, 'file') ~= 2
        error('Missing SPM CSD LFP data fixture: %s', data_file);
    end

    fit_dir = fullfile(work_root, [test_id '_' model]);
    if exist(fit_dir, 'dir'), rmdir(fit_dir, 's'); end
    mkdir(fit_dir);
    result.output_dir = fit_dir;

    D = spm_eeg_load(data_file);
    Dcopy = copy(D, fullfile(fit_dir, 'LFP_BG.mat'));

    DCM = dts_meeg_csd_lfp_template(model, Dcopy.fullfile(), cfg);
    DCM.name = fullfile(fit_dir, ['DCM_' test_id '_' model '_LFP']);
    DCM = spm_dcm_csd(DCM);

    residual = spm_vec(DCM.Rc);
    observed = spm_vec(DCM.Hc) + residual;
    observed_power = sum(abs(observed).^2);
    residual_power = sum(abs(residual).^2);

    result.F = double(full(DCM.F));
    result.R2 = 100*(observed_power - residual_power)/observed_power;
    result.Ep = spm_vec(DCM.Ep);
    result.file = [DCM.name '.mat'];
    result.ok = isfinite(result.F) && isfinite(result.R2) && all(isfinite(result.Ep));
    clear cleanup_threads
catch ME
    result.error = getReport(ME, 'extended', 'hyperlinks', 'off');
end
end
