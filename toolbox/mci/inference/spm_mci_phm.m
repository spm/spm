function [logev] = spm_mci_phm (L)
% Compute Log Evidence using Posterior Harmonic Mean (PHM)
% FORMAT [logev] = spm_mci_phm (L)
%
% L          [S x 1] vector containing log-likelihood of samples
% logev      log evidence from PHM
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

Lbar=mean(L);
logev=-log(mean(exp(Lbar-L)))+Lbar;