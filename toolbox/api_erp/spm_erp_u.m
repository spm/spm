function [U] = spm_erp_u(P,t)
% returns the [scalar] input for EEG models
% FORMAT [U] = spm_erp_u(P,t)
%
% P      - parameter structure
%   P.R  - scaling of [Poisson] parameter
%   P.N  - [DCT] parameter
%   P.U  - trial duration (seconds)
%
% t      - PST (seconds)
%
% See spm_fx_erp.m and spm_erp_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

global M

% get onset
%--------------------------------------------------------------------------
try
    onset = M.onset;
catch
    onset = 60;                                       % default of 60ms
end

% stimulus - subcortical impluse
%--------------------------------------------------------------------------
R  = exp(P.R).*[onset 1];                             % [Gamma] parameters
U  = spm_Gpdf(t*1000,R(1)*R(2),R(2))*1000;            % input

% Endogenous fluctuations (if P.N is specified)
%--------------------------------------------------------------------------
try
   f  = P.N(3);
   U  = U + P.N(1)*cos(2*pi*f*t) + P.N(2)*sin(2*pi*f*t);
end

