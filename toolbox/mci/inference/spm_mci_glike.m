function [L,e,st] = spm_mci_glike (P,M,U,Y,G)
% Gaussian Log-likelihood 
% FORMAT [L,e,st] = spm_mci_glike (P,M,U,Y,G)
%
% P         Parameters
% M         Model structure
% U         Inputs
% Y         Data
% G         Predictions (computed if not provided)
% 
% L         Log Likelihood
% e         Errors
% st        Status flag (0 for OK, -1 for problem)
%
% Assumes diagonal error covariance M.Ce
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_glike.m 6275 2014-12-01 08:41:18Z will $

st=0;
try
    G = G;
catch
    [G,tmp,st] = spm_mci_fwd (P,M,U);
end

e = Y-G;
L = -0.5*trace(M.iCe*e'*e) + M.logdet_Ce - 0.5*M.N*log(2*pi);

if isnan(L)
    L=realmin;
end
