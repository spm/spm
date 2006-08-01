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
    O = M.onset;
catch
    O = 120/8;                                     % defaolt prior of 120ms
end

% stimulus parameters
%--------------------------------------------------------------------------
R     = exp(P.R).*[O 1];                           % [Gamma] parameters
U     = spm_Gpdf(t*1000,R(1)*R(2),R(2))*1000;      % input
for i = 1:length(P.N)
    U = U + P.N(i)*cos((i - 1)*pi*t/P.U);          % [DCT] fluctuations
end
