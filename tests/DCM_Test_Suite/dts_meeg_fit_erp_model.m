function result = dts_meeg_fit_erp_model(test_id, model, data_file, spatial, generator_model, work_root, cfg)
% Fit one ERP-family DCM model to one staged M/EEG dataset.
% Authored by Pranay Yadav in 2026

estimation_mode = 'unspecified';
if isfield(cfg, 'estimation_mode')
    estimation_mode = cfg.estimation_mode;
end

result = struct( ...
    'test_id', test_id, ...
    'model', model, ...
    'spatial', spatial, ...
    'generator_model', generator_model, ...
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

    if exist(data_file, 'file') ~= 2
        error('Missing M/EEG input fixture: %s', data_file);
    end

    fit_dir = fullfile(work_root, [test_id '_' model]);
    if exist(fit_dir, 'dir'), rmdir(fit_dir, 's'); end
    mkdir(fit_dir);
    result.output_dir = fit_dir;

    D = spm_eeg_load(data_file);
    [~, data_name, data_ext] = fileparts(data_file);
    Dcopy = copy(D, fullfile(fit_dir, [data_name data_ext]));

    try
        gain = D.inv{D.val}.gainmat;
        if ~isempty(gain)
            [pth, nam, ext] = fileparts(gain);
            if isempty(pth), pth = D.path; end
            src = fullfile(pth, [nam ext]);
            if exist(src, 'file')
                copyfile(src, fullfile(fit_dir, [nam ext]));
            end
        end
    catch
    end

    switch spatial
        case 'LFP'
            DCM = dts_meeg_lfp_template(model, Dcopy.fullfile(), cfg);
            if strcmp(test_id, 'test_05_sample_lfp')
                DCM.xU.X = [0; 1];
                DCM.xU.name = {'deviant'};
                DCM.options.trials = [1 2];
            end
        case {'IMG','ECD'}
            DCM = dts_meeg_sensor_template(model, Dcopy.fullfile(), spatial, cfg);
            if strcmp(spatial, 'IMG')
                DCM.M.dipfit.Nm = cfg.IMG_Nm;
                DCM.M.dipfit.radius = cfg.IMG_radius;
            else
                DCM.options.location = cfg.ECD_location;
            end
        otherwise
            error('Unknown spatial model: %s', spatial);
    end

    DCM.testing.test_id = test_id;
    DCM.testing.generator_model = generator_model;
    DCM.name = fullfile(fit_dir, ['DCM_' test_id '_' model '_' spatial]);
    DCM = spm_dcm_erp(DCM);

    predicted_power = 0;
    residual_power = 0;
    for t = 1:numel(DCM.R)
        residual = full(DCM.R{t});
        fitted = full(DCM.H{t});
        predicted_power = predicted_power + sum(abs(fitted(:)).^2);
        residual_power = residual_power + sum(abs(residual(:)).^2);
    end

    result.F = double(full(DCM.F));
    result.R2 = 100*predicted_power/(predicted_power + residual_power);
    result.Ep = spm_vec(DCM.Ep);
    result.file = [DCM.name '.mat'];
    result.ok = isfinite(result.F) && isfinite(result.R2) && all(isfinite(result.Ep));
    clear cleanup_threads
catch ME
    result.error = getReport(ME, 'extended', 'hyperlinks', 'off');
end
end
