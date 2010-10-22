function [m,v] = spm_csd_chf(P,M,U)
% Characteristic (expected) frequency of a NMM
% FORMAT [G,w] = spm_csd_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% m - expected frequency
% v - dispersion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_chf.m 4095 2010-10-22 19:37:51Z karl $
 
 
% compute log-spectral density (without noise)
%==========================================================================
P.a(1) = -0;
P.a(2) = -32;
P.b    = P.b - 32;
P.c    = P.c - 32;
 
[G,w]  = spm_csd_mtf(P,M,U);
 
% compute moments (treating abs(G{i}) as a density)
%--------------------------------------------------------------------------
for i = 1:length(G)
    p     = abs(G{i});
    p     = p/sum(p);
    m(i)  = w'*p;
    v(i)  = ((w - m(i)).^2)'*p;
end
 
% average over trials
%--------------------------------------------------------------------------
m    = mean(m);
v    = mean(v);
