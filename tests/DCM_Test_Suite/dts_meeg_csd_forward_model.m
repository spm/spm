function result = dts_meeg_csd_forward_model(model)
% Run prior CSD forward passes for one- and four-node LFP networks.
% Authored by Pranay Yadav in 2026

result = struct( ...
    'test_id', 'test_06_csd_forward', ...
    'model', model, ...
    'spatial', 'LFP', ...
    'output_norm', NaN, ...
    'architecture_norms', [], ...
    'ok', false, ...
    'error', '');

try
    architectures = {'single','four'};
    output_norms = NaN(size(architectures));

    for a = 1:numel(architectures)
        architecture = architectures{a};

        switch architecture
            case 'single'
                A = {0, 0, 0};
                B = {0};
                C = 1;
                X = sparse(1, 0);
            case 'four'
                A = {zeros(4), zeros(4), zeros(4)};
                A{1}(2, 1) = 1;
                A{1}(4, 3) = 1;
                A{2}(1, 2) = 1;
                A{2}(3, 4) = 1;
                A{3}(3, 1) = 1;
                A{3}(1, 3) = 1;
                A{3}(4, 2) = 1;
                A{3}(2, 4) = 1;
                B = {double((A{1} + A{2} + A{3}) > 0)};
                C = [1; 0; 1; 0];
                X = sparse([0; 1]);
            otherwise
                error('Unknown CSD forward architecture: %s', architecture);
        end

        [P, pC] = spm_dcm_neural_priors(A, B, C, model);
        Ns = size(A{1}, 1);
        dipfit = struct('model', model, 'type', 'LFP', ...
            'location', 0, 'Ns', Ns, 'Nc', Ns);
        [P, pC] = spm_L_priors(dipfit, P, pC);
        P = spm_ssr_priors(P, pC);
        [x, f] = spm_dcm_x_neural(P, model);

        M = [];
        M.model = model;
        M.f = f;
        M.g = 'spm_gx_erp';
        M.x = x;
        M.n = length(spm_vec(x));
        M.m = Ns;
        M.u = sparse(Ns, 1);
        M.l = Ns;
        M.Hz = (2:64)';
        M.dt = 1/128;
        M.dipfit = dipfit;

        U = [];
        U.X = X;

        y = spm_csd_mtf(P, M, U);
        y_vec = spm_vec(y);
        output_norms(a) = double(full(norm(y_vec)));

        if isempty(y_vec) || ~all(isfinite(y_vec)) || output_norms(a) <= eps
            error('%s CSD output is empty, non-finite, or near-zero.', architecture);
        end
    end

    result.output_norm = min(output_norms);
    result.architecture_norms = output_norms;
    result.ok = true;
catch ME
    result.error = getReport(ME, 'extended', 'hyperlinks', 'off');
end
end
