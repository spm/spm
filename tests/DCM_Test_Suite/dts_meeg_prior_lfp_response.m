function y = dts_meeg_prior_lfp_response(model, dt, ns, integrator)
% Generate the default-prior one-node LFP response used by simulated tests.
% Authored by Pranay Yadav in 2026

if nargin < 4 || isempty(integrator)
    integrator = 'spm_int_L';
end

A = {0, 0, 0};
B = {0};
C = 1;
pE = spm_dcm_neural_priors(A, B, C, model);
[x, f] = spm_dcm_x_neural(pE, model);

M = [];
M.model = model;
M.f     = f;
M.x     = x;
M.n     = length(spm_vec(x));
M.m     = size(pE.C, 2);
M.ns    = ns;
M.ons   = 60;
M.dur   = 16;
M.integrator = integrator;

U = [];
U.dt = dt;
U.X  = [];

states = spm_gen_erp(pE, M, U);
dipfit = struct('model', model, 'type', 'LFP', 'location', 0, ...
    'Ns', 1, 'Nc', 1);
gE = spm_L_priors(dipfit);
L  = spm_lx_erp(gE, dipfit);
x0 = ones(size(states{1}, 1), 1)*spm_vec(M.x)';
y  = full((states{1} - x0)*L');
end
