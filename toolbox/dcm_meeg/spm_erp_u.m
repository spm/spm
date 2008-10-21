function [u] = spm_erp_u(t,P,M)
% returns the [scalar] input for EEG models (Gaussian function)
% FORMAT [u] = spm_erp_u(t,P,M)
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
% $Id: spm_erp_u.m 2374 2008-10-21 18:52:29Z karl $


% stimulus - subcortical impulse
%--------------------------------------------------------------------------
nu    = length(M.ons);
u     = sparse(length(t),nu);
t     = t*1000;
for i = 1:nu
   delay  = M.ons(i) + P.R(i,1)*128;
   scale  = exp(P.R(i,2))*64;
   u(:,i) = exp(-(t - delay).^2/scale)*32;
end

