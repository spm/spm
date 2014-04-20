function [u] = spm_erp_u(t,P,M)
% returns the [scalar] input for EEG models (Gaussian function)
% FORMAT [u] = spm_erp_u(t,P,M)
% t      - PST (seconds)
% P      - parameter structure
%   P.R  - scaling of [Gaussian] parameters
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
% $Id: spm_erp_u.m 5964 2014-04-20 09:48:58Z karl $


% preliminaries - check durations (ms)
%--------------------------------------------------------------------------
try
    if length(M.dur) ~= length(M.ons)
        M.dur = M.dur(1) + M.ons - M.ons;
    end
catch
    M.dur = 32 + M.ons - M.ons;
end

% check sustained input (0,1)
%--------------------------------------------------------------------------
try
    if length(M.sus) ~= length(M.ons)
        M.sus = M.sus(1) + M.ons - M.ons;
    end
catch
    M.sus = 0 + M.ons - M.ons;
end

% stimulus – Gaussian (subcortical) impulse
%--------------------------------------------------------------------------
nu    = length(M.ons);
u     = sparse(length(t),nu);
t     = t*1000;
for i = 1:nu
    
    % Gaussian bump function
    %----------------------------------------------------------------------
    delay  = M.ons(i) + 128*P.R(i,1);
    scale  = M.dur(i) * exp(P.R(i,2));
    U      = exp(-(t - delay).^2/(2*scale^2));
    
    % sustained inputs
    %----------------------------------------------------------------------
    try
        prop = M.sus(i)*exp(P.R(i,3));
    catch
        prop = M.sus(i);
    end
    U      = prop*cumsum(U)/sum(U) + U*(1 - prop);
    u(:,i) = 32*U;
end



