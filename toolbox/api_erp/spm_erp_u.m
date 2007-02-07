function [U,N] = spm_erp_u(P,t)
% returns the [scalar] input for EEG models
% FORMAT [U,N] = spm_erp_u(P,t)
%
% P      - parameter structure
%   P.R  - scaling of [Gamma] parameters
%   P.N  - [DCT] parameter[s]
%
% t      - PST (seconds)
%
% U   - stimulus-related (subcortical) input
% N   - non-specifc fluctuations
%
% See spm_fx_erp.m and spm_erp_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% stimulus - subcortical impulse: [Gamma] parameters
%--------------------------------------------------------------------------
U     = 0;
for i = 1:length(P.ons)
   magni = P.M(i);
   shape = P.ons(i)*exp(P.R(i,1));
   scale = exp(P.R(i,2));
   U     = U + magni*spm_Gpdf(t*1000,shape,scale)*1000;
end

% Endogenous fluctuations (if P.N is specified)
%--------------------------------------------------------------------------
try
   for i = 1:size(P.N,1)
      N(i,:) = P.N(i,1)*cos(2*pi*P.N(i,3)*t) + P.N(i,2)*sin(2*pi*P.N(i,3)*t);
   end
end

