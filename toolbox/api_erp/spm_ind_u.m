function [U,B] = spm_ind_u(t,P,M)
% returns the [scalar] input for EEG models
% FORMAT [U,B] = spm_ind_u(t,P,M)
%
% P      - parameter structure
%   P.R  - scaling of [Gamma] parameters
%   P.N  - [DCT] parameter[s]
%
% t      - PST (seconds)
%
% U   - stimulus-related (subcortical) input
% B   - non-specifc background fluctuations
%
% See spm_fx_ind.m and spm_ind_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% stimulus - subcortical impulse
%--------------------------------------------------------------------------
ampli = exp(P.R(1))*32;
delay = M.ons(1)*exp(P.R(2));
scale = exp(P.R(3))*32;

U     = ampli*exp(-(t*1000 - delay).^2/scale);

% Endogenous fluctuations (if P.N is specified)
%--------------------------------------------------------------------------
w      = exp(j*6.2832*P.N(3)*t);
B      = P.N(1:2)*[real(w); imag(w)];


