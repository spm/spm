function [u] = spm_ind_u(t,P,M)
% returns the [scalar] input for EEG models (Gamma function)
% FORMAT [u] = spm_ind_u(t,P,M)
%
% P     - parameter structure
%   P.R - input parameters
%
% t     - PST (seconds)
%
% u   - stimulus-related (subcortical) input
%
% See spm_erp_u, spm_fx_ind.m and spm_ind_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ind_u.m 1186 2008-03-05 12:52:57Z karl $

% stimulus - subcortical impulse
%--------------------------------------------------------------------------
nu    = length(M.ons);
u     = sparse(length(t),nu);
t     = t*1000;

% Gaussian form
%--------------------------------------------------------------------------
for i = 1:nu
   delay  = M.ons(i) + P.R(i,1)*128;
   scale  = exp(P.R(i,2))*2048;
   u(:,i) = exp(-(t - delay).^2/scale);
end

return

% NB: Code for Gamma function input
%==========================================================================

% stimulus - subcortical impulse - a gamma function
%--------------------------------------------------------------------------
for i = 1:nu
    m      = P.R(i,1) + M.ons(i);
    v      = exp(P.R(i,2))*1024;
    u(:,i) = Gpdf(t,m*m/v,m/v);
end


function f = Gpdf(x,a,b)
%--------------------------------------------------------------------------
Q    = find(x > 0);
f    = x*0;
f(Q) = x(Q).^(a - 1).*exp(-b*x(Q))*(b^a)/gamma(a);


% stimulus - subcortical impulse
%--------------------------------------------------------------------------
nu    = length(M.ons);
u     = sparse(length(t),nu);
t     = t*1000;
for i = 1:nu
   delay  = exp(P.R(i,1))*M.ons(i);
   scale  = exp(P.R(i,2))*64;
   u(:,i) = exp(-(t - delay).^2/scale)*32;
end
