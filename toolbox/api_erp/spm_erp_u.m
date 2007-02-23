function [U,B] = spm_erp_u(t,P,M)
% returns the [scalar] input for EEG models
% FORMAT [U,B] = spm_erp_u(t,P,M)
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
% See spm_fx_erp.m and spm_erp_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% stimulus - subcortical impulse
%--------------------------------------------------------------------------
U     = sparse(1,length(t));
tms   = t*1000;
for i = 1:length(M.ons)
   delay = M.ons(i)*exp(P.R(i,2));
   scale = exp(P.R(i,3))*64;
   U     = U + 32*P.R(i,1)*exp(-(tms - delay).^2/scale);
end

% Endogenous fluctuations (if P.N is specified)
%--------------------------------------------------------------------------
n     = size(P.N,1);
B     = sparse(n,length(t));
try
    for i = 1:n
       w      = exp(j*6.2832*P.N(i,3)*t);
       B(i,:) = P.N(i,1:2)*[real(w); imag(w)];
    end
end

