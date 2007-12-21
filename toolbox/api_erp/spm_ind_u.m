function [U] = spm_ind_u(t,P,M)
% returns the [scalar] input for EEG models
% FORMAT [U] = spm_ind_u(t,P,M)
%
% P     - parameter structure
%   P.R - input parameters
%
% t     - PST (seconds)
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
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id$

% stimulus - subcortical impulse - a gamma function
%--------------------------------------------------------------------------
ons   = M.ons(:);
dur   = M.dur*1000;
t     = t(:)';
j     = 1:length(t);
for i = 1:length(ons)
    m      = exp(P.R(i,1))* ons(i);
    v      = exp(P.R(i,2))*(dur/16)^2;
    U(i,:) = Gpdf(t*1000,m*m/v,m/v);
end

function f = Gpdf(x,h,l)
%--------------------------------------------------------------------------
Q     = find(x > 0);
f     = x*0;
f(Q)  = x(Q).^(h - 1).*exp(-l*x(Q))*(l^h)/gamma(h);