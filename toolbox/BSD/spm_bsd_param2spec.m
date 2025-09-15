function spec = spm_bsd_param2spec(P, M)
% Convert BSD model parameters to spectral parameters structure
% FORMAT spec = spm_bsd_param2spec(P, M)
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

% Johan Medrano

% Initialize output structure
spec = []; 

% Map log-amplitude parameters to amplitude
if isfield(P, 'a')
    spec.ampl = real(exp(P.a));
end

% Map frequency parameters to frequency within prior bounds
if isfield(P, 'f')
    spec.freq = real(M.pV.f(:, 1) + (0.5 + 0.5 * tanh(P.f)) ... 
        .* (M.pV.f(:, 2) - M.pV.f(:, 1)));
end

% Map log-FWHM parameters to FWHM using prior scaling
if isfield(P, 'S')
    spec.fwhm = real(sqrt(exp(P.S) .* M.pV.S * 8 * log(2))); 
end

% Map log-noise parameters to noise using prior scaling
if isfield(P, 'b')
    spec.noise = real(exp(P.b) .* M.pV.b);
end


end