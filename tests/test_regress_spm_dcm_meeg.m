classdef test_regress_spm_dcm_meeg < matlab.unittest.TestCase
% SPM DCM M/EEG regression tests
% Authored by Pranay Yadav in 2026

% Get options for running tests
properties (Constant)
    Options = test_regress_spm_dcm_meeg_options();
end

% Test parameterisation based on model, integration scheme and data cases
properties (TestParameter)

	% Each test will always loop over all models
    Model = {'ERP','SEP','LFP','CMC','CMM','NMM','MFM','NMDA','CMM_NMDA'};
    
    % Forward simulations will be tested for integrators listed here
    Integrator = {'spm_int_L','spm_int_B','spm_int_J'};
    
    % DCM-ERP inversions are tested against multiple data cases listed here
    FixedDataCase = {'simulated_evoked_lfp','synthetic_triphasic_lfp', ...
        'sample_lfp','sample_sensor_img','sample_sensor_ecd'};
end

% Base prep: paths and env
methods (TestClassSetup)
    function setupSPM(~)
    
    	% Run in single thread only (for max determinism)
        test_regress_spm_dcm_meeg_thread_state(maxNumCompThreads());
        maxNumCompThreads(1);
        test_dir = fileparts(mfilename('fullpath'));
        support_dir = fullfile(test_dir, 'DCM_Test_Suite');
        addpath(support_dir);
        dts_meeg_setup_paths();
        
        % Configure settings based on options specified above
        cfg = dts_meeg_test_defaults(test_regress_spm_dcm_meeg.Options);
        
        % Generate test fixtures with fwd model for sensor-level ECD & IMG fits
        dts_meeg_sample_sensor_fixture(cfg, 'testrunner');
    end
end

% Restore num threads to default once all tests are done
methods (TestClassTeardown)
    function teardownSPM(~)
        original_num_threads = test_regress_spm_dcm_meeg_thread_state();
        if ~isempty(original_num_threads)
            maxNumCompThreads(original_num_threads);
        end
    end
end

% Test for model stability by examining eigvals of steady-state Jacobian
methods (Test, TestTags = {'stability'})
    function test_regress_01_steady_state_jacobian(testCase, Model)
        cfg = testCase.config_for_model(Model);
        [P, M] = dts_meeg_neural_fixture(Model, 'single', 'spm_int_L', 1000, 0.001);
        u = sparse(M.m, 1);
        f0 = full(spm_vec(feval(M.f, M.x, u, P, M)));
        J = full(spm_cat(spm_diff(M.f, M.x, u, P, M, 1)));
        max_real_eig = max(real(eig(J)));

        testCase.verifyTrue(isfinite(norm(f0)), ...
            sprintf('%s non-finite gradient norm', Model));
        testCase.verifyLessThan(max_real_eig, cfg.thresholds.jacobian_max_real_eig, ...
            sprintf('%s max real eig %.15g >= %.15g', ...
            Model, max_real_eig, cfg.thresholds.jacobian_max_real_eig));
    end
end

% Test for successful initialisation of model to steady state solution
methods (Test, TestTags = {'init'})
    function test_regress_02_initialisation_baseline(testCase, Model)
        cfg = testCase.config_for_model(Model);
        [P, M, U] = dts_meeg_neural_fixture(Model, 'single', 'spm_int_L', 1000, 0.001);
        u0 = sparse(M.m, 1);
        f0 = full(spm_vec(feval(M.f, M.x, u0, P, M))); 
        y = spm_gen_erp(P, M, U); % simulate ERP, this does initialisation
        states = full(y{1});
        x0 = spm_vec(M.x)';
        baseline = ones(size(states, 1), 1)*x0;
        centered = states - baseline;
        scale = max(1, max(abs(centered(:))));
        start_error = norm(centered(1, :))/scale;
        end_error = norm(centered(end, :))/scale;

        testCase.verifyTrue(isfinite(norm(f0)) && isfinite(start_error) && isfinite(end_error), ...
            sprintf('%s non-finite baseline diagnostic', Model));
        testCase.verifyLessThan(start_error, cfg.thresholds.baseline_error, ...
            sprintf('%s start baseline %.15g >= %.15g', ...
            Model, start_error, cfg.thresholds.baseline_error));
        testCase.verifyLessThan(end_error, cfg.thresholds.baseline_error, ...
            sprintf('%s end baseline %.15g >= %.15g', ...
            Model, end_error, cfg.thresholds.baseline_error));
    end
end

% Test if all models can be integrated successfully for 1-node & 4-node network
methods (Test, TestTags = {'integrate'})
    function test_regress_03_integrator_forward(testCase, Model, Integrator)
        cfg = testCase.config_for_model(Model);
        architectures = {'single','four'};
		
		% Each integrator is compared against the *J integrator
        for a = 1:numel(architectures)
            arch = architectures{a};
            [P, M, U] = dts_meeg_neural_fixture(Model, arch, 'spm_int_J', 1000, 0.001);
            ref = double(full(spm_vec(spm_gen_erp(P, M, U))));

            [P, M, U] = dts_meeg_neural_fixture(Model, arch, Integrator, 1000, 0.001);
            actual = double(full(spm_vec(spm_gen_erp(P, M, U))));
            c = corrcoef(ref(:), actual(:));
            r = c(1, 2);
            testCase.verifyGreaterThanOrEqual(r, cfg.thresholds.integrator_min_corr, ...
                sprintf('%s %s %s corr %.15g < %.15g', ...
                Model, arch, Integrator, r, cfg.thresholds.integrator_min_corr));
        end
    end
end

% Test if each model can fit ERP generated by its own default priors
methods (Test, TestTags = {'self_erp_fit'})
    function test_regress_04_self_generated_lfp_fit(testCase, Model)
        test_id = 'test_04_self_generated_lfp';
        cfg = testCase.config_for_fit(Model, '');

        if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
        data_dir = fullfile(cfg.output_root, [test_id '_' Model '_generated_lfp']);
        if exist(data_dir, 'dir'), rmdir(data_dir, 's'); end
        mkdir(data_dir);

        data_file = fullfile(data_dir, ['self_generated_lfp_' Model '.mat']);
        rng_state = rng;
        cleanup = onCleanup(@()rng(rng_state));
        y = dts_meeg_prior_lfp_response(Model, 0.001, 500, cfg.integrator);
        y = y(:) - y(1);
        y = y/std(y);
        rng(1, 'twister'); % Fix random noise 
        y = y + randn(size(y))*exp(-6/2);
        D = dts_meeg_write_lfp_fixture(y, 0.001, data_file, ['simulated_' Model]);
        data_file = D.fullfile();
        spm_eeg_load(data_file);
        clear cleanup

        result = dts_meeg_fit_erp_model(test_id, Model, data_file, 'LFP', Model, cfg.output_root, cfg);
        testCase.verify_fit_completed(result);
        testCase.verifyGreaterThanOrEqual(result.R2, cfg.thresholds.self_lfp_min_r2, ...
            sprintf('%s %s R2 %.15g < %.15g', ...
            test_id, Model, result.R2, cfg.thresholds.self_lfp_min_r2));

        fit_passed = isempty(result.error) && result.ok && result.R2 >= cfg.thresholds.self_lfp_min_r2;
        if ~cfg.keep_outputs && (fit_passed || ~cfg.preserve_on_failure)
            if isfield(result, 'output_dir') && exist(result.output_dir, 'dir')
                rmdir(result.output_dir, 's');
            end
            if exist(data_dir, 'dir')
                rmdir(data_dir, 's');
            end
        end
    end
end

% Test if each model can fit to various data cases
methods (Test, TestTags = {'data_erp_fit'})
    function test_regress_05_fixed_data_fit(testCase, Model, FixedDataCase)
        data_case = dts_meeg_fixed_data_cases(FixedDataCase);
        cfg = testCase.config_for_fit(Model, data_case.name);
        data_file = fullfile(cfg.prep_root, data_case.prep_dir, data_case.file_name);

        testCase.assertEqual(exist(data_file, 'file'), 2, ...
            sprintf(['Missing prepared DCM M/EEG fixture:\n%s\n' ...
            'Run tests/DCM_Test_Suite/dts_prep01_make_test_fixtures.m.'], data_file));

        if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
        result = dts_meeg_fit_erp_model(data_case.test_id, Model, data_file, ...
            data_case.spatial, data_case.generator_model, cfg.output_root, cfg);
        testCase.verify_fit_completed(result);

        errors = dts_meeg_check_fit_reference(result, cfg);
        fit_passed = isempty(result.error) && result.ok && isempty(errors);
        if ~cfg.keep_outputs && (fit_passed || ~cfg.preserve_on_failure)
            if isfield(result, 'output_dir') && exist(result.output_dir, 'dir')
                rmdir(result.output_dir, 's');
            end
        end
        if ~isempty(errors)
            testCase.verifyFail(sprintf('%s %s reference check failed:\n%s', ...
                data_case.test_id, Model, strjoin(errors, newline)));
        end
    end
end

% Test if CSD forward simulations work for 1-node and 4-node networks
methods (Test, TestTags = {'csd_gen'})
    function test_regress_06_csd_forward(testCase, Model)
        result = dts_meeg_csd_forward_model(Model);

        if ~isempty(result.error)
            testCase.verifyFail(sprintf('%s %s error:\n%s', ...
                result.test_id, result.model, result.error));
            return
        end

        testCase.verifyTrue(result.ok, sprintf('%s %s did not produce finite CSD output', ...
            result.test_id, result.model));
        testCase.verifyGreaterThan(result.output_norm, eps, ...
            sprintf('%s %s produced near-zero CSD output', result.test_id, result.model));
    end
end

% Test if models can be fit to CSD from sample LFP data
methods (Test, TestTags = {'data_csd_fit'})
    function test_regress_07_csd_lfp_fit(testCase, Model)
        cfg = testCase.config_for_fit(Model, '');
        if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
        result = dts_meeg_fit_csd_model(Model, cfg.output_root, cfg);
        testCase.verify_fit_completed(result);

        errors = dts_meeg_check_fit_reference(result, cfg);
        fit_passed = isempty(result.error) && result.ok && isempty(errors);
        if ~cfg.keep_outputs && (fit_passed || ~cfg.preserve_on_failure)
            if isfield(result, 'output_dir') && exist(result.output_dir, 'dir')
                rmdir(result.output_dir, 's');
            end
        end
        if ~isempty(errors)
            testCase.verifyFail(sprintf('%s %s reference check failed:\n%s', ...
                result.test_id, Model, strjoin(errors, newline)));
        end
    end
end

methods
	% Get default options and apply custom ones for test
    function cfg = config_for_model(~, model)
        options = test_regress_spm_dcm_meeg.Options;
        options.models = {model};
        cfg = dts_meeg_test_defaults(options);
    end
	
	% Get customized options for test and adjust Nmax based on fixture
    function cfg = config_for_fit(testCase, model, fixed_data_case)
        cfg = testCase.config_for_model(model);
        if any(strcmp(cfg.strict_fixed_data_cases, fixed_data_case))
            cfg.estimation_mode = 'strict';
            cfg.Nmax = cfg.strict_fit_Nmax;
            cfg.thresholds.min_R2 = cfg.thresholds.strict_min_R2;
        else
            cfg.estimation_mode = 'basic';
            cfg.Nmax = cfg.basic_fit_Nmax;
            cfg.thresholds.min_R2 = cfg.thresholds.basic_min_R2;
        end
    end
	
	% Check if fitting tests ran
    function verify_fit_completed(testCase, result)
        if ~isempty(result.error)
            testCase.verifyFail(sprintf('%s %s error:\n%s', ...
                result.test_id, result.model, result.error));
            return
        end

        testCase.verifyTrue(result.ok, sprintf('%s %s did not produce finite F/R2/Ep', ...
            result.test_id, result.model));
        testCase.verifyTrue(isfinite(result.F), sprintf('%s %s F is not finite', ...
            result.test_id, result.model));
        testCase.verifyTrue(isfinite(result.R2), sprintf('%s %s R2 is not finite', ...
            result.test_id, result.model));
        testCase.verifyTrue(all(isfinite(result.Ep)), sprintf('%s %s Ep is not finite', ...
            result.test_id, result.model));
    end
end
end

% Options for this test
function options = test_regress_spm_dcm_meeg_options()
options = [];

% Whether the test runner is expected to execute tests in parallel
options.parallel = false;  % Only used if running tests manually

% Keep DCM output files after successful test execution
options.keep_outputs = false;

% Preserve generated outputs when a fit fails
options.preserve_on_failure = true;

% Use deterministic output paths under the test-data tree
options.deterministic_work_root = true;

% NLSI iterations for fixed data cases that need strict estimation ("strict")
options.strict_fit_Nmax = 32;

% NLSI iterations for broad regression fits that only check inversion path ("basic")
options.basic_fit_Nmax = 4;

% Names of FixedDataCase tests that run in "strict" mode for all models
options.strict_fixed_data_cases = {'simulated_evoked_lfp', 'sample_lfp'};

% Stability bound: largest real Jacobian eigval must be less than this
options.thresholds.jacobian_max_real_eig = 1e-9;

% Maximum normalised deviation from a model's initialised baseline states
options.thresholds.baseline_error = 0.1;

% Minimum correlation of model output with the spm_int_J reference output
options.thresholds.integrator_min_corr = 0.75;

% Optional R2 floor for self-gen LFP fits (0 disables this)
options.thresholds.self_lfp_min_r2 = 0;

% Allowed absolute drift in free energy 
options.thresholds.F_margin = 5;

% Allowed maximum absolute drift in posterior estimates
options.thresholds.Ep_margin = 0.1;

% Allowed absolute R2% drift from reference fits
options.thresholds.R2_margin = 5;

% Optional absolute R2% floor for "strict" fits (0 disables it)
options.thresholds.strict_min_R2 = 0;

% Optional absolute R2% floor for "basic" fits (0 disables it)
options.thresholds.basic_min_R2 = 0;
end

function value = test_regress_spm_dcm_meeg_thread_state(value)
persistent original_num_threads
if nargin > 0
    original_num_threads = value;
end
value = original_num_threads;
end
