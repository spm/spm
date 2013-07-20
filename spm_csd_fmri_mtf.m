function [y,w,s] = spm_csd_fmri_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_fmri_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects (induces expansion around steady state)
%
% y - y(N,nc,nc} - cross-spectral density for nc channels {trials}
%                - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - directed transfer functions (complex)
%
% When called with U this function will return a cross-spectral response
% for each of the condition-specific parameters specified in U.X; otherwise
% it returns the complex CSD for the parameters in P.
%
% When the observer function M.g is specified the CSD response is
% supplemented with channel noise in sensor space; otherwise the CSD
% pertains to hidden states.
%
% NB: requires M.u to specify the number of endogenous inputs
% This routine and will solve for the (hidden) steady state and use it as
% the expansion point for subsequent linear systems analysis (if trial
% specific effects are specified).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_fmri_mtf.m 5587 2013-07-20 15:37:17Z karl $
 

% compute log-spectral density
%==========================================================================
 
% frequencies of interest
%--------------------------------------------------------------------------
w    = M.Hz;
nw   = length(w);
 
% number of regions and exogenous (neuronal) inputs or sources
%--------------------------------------------------------------------------
nc   = M.l;
ns   = length(M.u);

% spectrum of neuronal innovations (Gu) and noise (Gs)
%==========================================================================
for i = 1:ns
    Gu(:,i) = exp(P.a(1,i))*w.^(-exp(P.a(2,i)));
end
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i) - 2)*w.^(-exp(P.c(2,i)));
end

% add exogenous input
%--------------------------------------------------------------------------
P.C   = eye(ns,ns);

 
% transfer functions (FFT of kernel)
%==========================================================================

S     = spm_dcm_mtf(P,M);

% [cross]-spectral density from neuronal innovations
%----------------------------------------------------------------------
G     = zeros(nw,nc,nc);
for i = 1:nc
    for j = 1:nc
        for k = 1:ns
            Gij      = S(:,i,k).*conj(S(:,j,k));
            G(:,i,j) = G(:,i,j) + Gij.*Gu(:,k);
        end
    end
end

% and add channel noise
%==========================================================================
if isfield(M,'g')
    for i = 1:nc
        G(:,i,i) = G(:,i,i) + Gs(:,i);
    end
    y = G + Gs;
else
    y = g;
end
