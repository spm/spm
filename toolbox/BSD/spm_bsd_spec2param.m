function P = spm_bsd_spec2param(spec, M, P)
% Convert spectral parameters structure to model parameters
% FORMAT P = spm_bsd_spec2param(spec, M, P)
%
% P   - structure containing model parameters:
%       .a : log-amplitude parameters
%       .f : frequency parameters
%       .S : log-FWHM (full width at half maximum) parameters
%       .b : log-noise parameters
%
% M   - structure containing prior values:
%       .pV.f : [Nx2] matrix, prior frequency bounds (Hz)
%       .pV.S : prior FWHM scaling
%       .pV.b : prior noise scaling
%
% spec - output structure containing spectral parameters:
%       .ampl : amplitude (exp-transformed from P.a)
%       .freq : frequency (transformed from P.f and prior bounds)
%       .fwhm : FWHM (derived from P.S and prior scaling)
%       .noise: noise (exp-transformed from P.b and prior scaling)
%
% This function maps the parameter estimates from a Bayesian Spectral
% Decomposition (BSD) model to a structure of spectral parameters for
% further analysis or model evaluation.
%
%_______________________________________________________________________
% Copyright (C) 2024-2025 Wellcome Trust Centre for Neuroimaging

try P; catch, P = []; end

% forward: a = exp(P.a);
if isfield(spec, 'ampl')
    if any(spec.ampl < 0), warning('Negative amplitude!'); end
    P.a = log(spec.ampl); 
end

% forward: f = M.pV.f(:, 1) + (0.5 + 0.5 * tanh(P.f)) ...
%               .* (M.pV.f(:, 2) - M.pV.f(:, 1));
if isfield(spec, 'freq')
    if any(spec.freq < 0), warning('Negative frequency!'); end
    P.f = atanh(2 * (spec.freq - M.pV.f(:, 1)) ... 
        ./ (M.pV.f(:, 2) - M.pV.f(:, 1)) - 1);
end

% forward: S = exp(P.S) .* M.pV.S; 
if isfield(spec, 'fwhm')
    if any(spec.fwhm < 0), warning('Negative FWHM!'); end
    S   = spec.fwhm.^2 ./ (8 * log(2));
    P.S = log(S) - log(M.pV.S);
end

% forward: b = exp(P.b) .* M.pV.b;
if isfield(spec, 'noise')
    P.b = log(abs(spec.noise)) - log(M.pV.b);
end
end

