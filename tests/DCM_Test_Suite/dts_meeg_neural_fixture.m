function [P, M, U] = dts_meeg_neural_fixture(model, architecture, integrator, ns, dt)
% Build a default-prior neural model for forward diagnostics.
% Authored by Pranay Yadav in 2026

if nargin < 2 || isempty(architecture), architecture = 'single'; end
if nargin < 3 || isempty(integrator), integrator = 'spm_int_L'; end
if nargin < 4 || isempty(ns), ns = 1000; end
if nargin < 5 || isempty(dt), dt = 0.001; end

switch lower(architecture)
    case {'single', 'one', 'one_node'}
        A = {0, 0, 0};
        B = {0};
        C = 1;
        X = sparse(1, 0);
    case {'four', 'four_node', 'four_nodes'}
        A = {zeros(4), zeros(4), zeros(4)};

        % Ipsilateral forward edges
        A{1}(2, 1) = 1;
        A{1}(4, 3) = 1;

        % Ipsilateral backward edges
        A{2}(1, 2) = 1;
        A{2}(3, 4) = 1;

        % Lateral homologous edges
        A{3}(3, 1) = 1;
        A{3}(1, 3) = 1;
        A{3}(4, 2) = 1;
        A{3}(2, 4) = 1;

        B = {double((A{1} + A{2} + A{3}) > 0)};
        C = [1; 0; 1; 0];
        X = sparse([0; 1]);
    otherwise
        error('Unknown diagnostic architecture: %s', architecture);
end

P = spm_dcm_neural_priors(A, B, C, model);
[x, f] = spm_dcm_x_neural(P, model);

M = [];
M.model = model;
M.f = f;
M.x = x;
M.m = size(P.C, 2);
M.n = length(spm_vec(x));
M.ons = 60;
M.dur = 16;
M.ns = ns;
M.integrator = integrator;

U = [];
U.dt = dt;
U.X = X;
end
