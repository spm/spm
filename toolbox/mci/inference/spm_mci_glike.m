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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_glike.m 6548 2015-09-11 12:39:47Z will $

st=0;
try
    G = G;
catch
    if isfield(M,'IS')
        [G,tmp,st] = feval(M.IS,P,M,U);
    else
        [G,tmp,st] = spm_mci_fwd (P,M,U);
    end
end

e = Y-G;
d = size(e,2);
L = -0.5*trace(M.iCe*e'*e) - 0.5*M.N*M.logdet_Ce - 0.5*M.N*d*log(2*pi);

if isnan(L)
    L=realmin;
end
