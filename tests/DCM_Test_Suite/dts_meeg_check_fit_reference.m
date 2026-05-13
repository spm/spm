function errors = dts_meeg_check_fit_reference(results, cfg)
% Check fitted DCM results against shipped deterministic references.
% Authored by Pranay Yadav in 2026

errors = {};

if exist(cfg.reference_file, 'file') ~= 2
    error(['Missing DCM M/EEG reference file:\n%s\n' ...
        'Run tests/DCM_Test_Suite/dts_prep02_generate_reference_data.m ' ...
        'and ship/copy refdata with the SPM test data.'], cfg.reference_file);
end
x = load(cfg.reference_file);
reference = x.reference;

for i = 1:numel(results)
    result = results(i);
    rows = reference.rows;
    idx = find(strcmp({rows.test_id}, result.test_id) & strcmp({rows.model}, result.model), 1);
    if isempty(idx)
        errors{end+1} = sprintf('No reference row found for %s/%s', result.test_id, result.model); %#ok<AGROW>
        continue
    end
    ref = rows(idx);
    result.F = double(full(result.F));
    result.R2 = double(full(result.R2));
    ref.F = double(full(ref.F));
    ref.R2 = double(full(ref.R2));
    result.Ep = double(full(result.Ep));
    ref.Ep = double(full(ref.Ep));

    if ~isempty(result.error)
        errors{end+1} = sprintf('%s %s error:\n%s', result.test_id, result.model, result.error); %#ok<AGROW>
        continue
    end
    if isfield(ref, 'Nmax') && isfield(result, 'Nmax') && double(full(result.Nmax)) ~= double(full(ref.Nmax))
        errors{end+1} = sprintf('%s %s Nmax mismatch: actual %.15g ref %.15g', ...
            result.test_id, result.model, double(full(result.Nmax)), double(full(ref.Nmax))); %#ok<AGROW>
    end
    F_margin = cfg.thresholds.F_margin;
    Ep_margin = cfg.thresholds.Ep_margin;
    R2_margin = cfg.thresholds.R2_margin;
    min_R2 = cfg.thresholds.min_R2;
    if isfield(ref, 'F_margin'), F_margin = double(full(ref.F_margin)); end
    if isfield(ref, 'Ep_margin'), Ep_margin = double(full(ref.Ep_margin)); end
    if isfield(ref, 'R2_abs_tol'), R2_margin = double(full(ref.R2_abs_tol)); end
    if isfield(ref, 'min_R2'), min_R2 = double(full(ref.min_R2)); end

    if ~isfinite(result.F)
        errors{end+1} = sprintf('%s %s F is not finite', result.test_id, result.model); %#ok<AGROW>
    elseif abs(result.F - ref.F) > F_margin
        errors{end+1} = sprintf('%s %s F drift: actual %.15g ref %.15g tol %.15g', ...
            result.test_id, result.model, result.F, ref.F, F_margin); %#ok<AGROW>
    end
    if ~isfinite(result.R2)
        errors{end+1} = sprintf('%s %s R2 is not finite', result.test_id, result.model); %#ok<AGROW>
    elseif min_R2 > 0 && (result.R2 < -eps || result.R2 > 100 + eps)
        errors{end+1} = sprintf('%s %s R2 out of percent range: %.15g', ...
            result.test_id, result.model, result.R2); %#ok<AGROW>
    else
        if abs(result.R2 - ref.R2) > R2_margin
            errors{end+1} = sprintf('%s %s R2 drift: actual %.15g ref %.15g tol %.15g', ...
                result.test_id, result.model, result.R2, ref.R2, R2_margin); %#ok<AGROW>
        end
        if min_R2 > 0 && result.R2 < min_R2
            errors{end+1} = sprintf('%s %s R2 below threshold: actual %.15g min %.15g', ...
                result.test_id, result.model, result.R2, min_R2); %#ok<AGROW>
        end
    end
    if numel(result.Ep) ~= numel(ref.Ep)
        errors{end+1} = sprintf('%s %s Ep length changed: actual %d ref %d', ...
            result.test_id, result.model, numel(result.Ep), numel(ref.Ep)); %#ok<AGROW>
    else
        ep_err = max(abs(result.Ep(:) - ref.Ep(:)));
        if ep_err > Ep_margin
            errors{end+1} = sprintf('%s %s Ep drift: max abs %.15g tol %.15g', ...
                result.test_id, result.model, ep_err, Ep_margin); %#ok<AGROW>
        end
    end
end
end
