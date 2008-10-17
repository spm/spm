function [u,dudP] = spm_erp_u_dudp(t,P,M)
% returns the [scalar] input for EEG models (Gaussian function)
% FORMAT [u,dudP] = spm_erp_u_dudp(t,P,M)
%
% P      - parameter structure
%   P.R  - scaling of [Gaussian] parameters
%
% t      - PST (seconds)
%
% u   - stimulus-related (subcortical) input
%
% See spm_fx_erp.m and spm_erp_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_erp_u_dudp.m 2353 2008-10-17 11:56:15Z karl $


% stimulus - subcortical impulse
%--------------------------------------------------------------------------
nu    = length(M.ons);
u     = sparse(length(t),nu);
t     = t*1000;
dudP = zeros(length(t),nu,nu*2);
for i = 1:nu
    delay  = M.ons(i) + P.R(i,1)*128;
    scale  = exp(P.R(i,2))*64;
    u(:,i) = exp(-(t - delay).^2/scale)*32;

    dudP(:,i,i) = 128*(t-M.ons(i)-P.R(i,1)*128)./exp(P.R(i,2))...
        .*exp(-1/64*(t-M.ons(i)-P.R(i,1)*128).^2./exp(P.R(i,2)));
    dudP(:,i,nu+i) = 1/2*(t-M.ons(i)-P.R(i,1)*128).^2./exp(P.R(i,2))...
        .*exp(-1/64*(t-M.ons(i)-P.R(i,1)*128).^2./exp(P.R(i,2)));

end

