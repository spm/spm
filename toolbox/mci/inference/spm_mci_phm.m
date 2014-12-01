function [logev] = spm_mci_phm (L)
% Compute Log Evidence using Posterior Harmonic Mean (PHM)
% FORMAT [logev] = spm_mci_phm (L)
%
% L          [S x 1] vector containing log-likelihood of samples
% logev      log evidence from PHM
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_phm.m 6275 2014-12-01 08:41:18Z will $

Lbar=mean(L);
logev=-log(mean(exp(Lbar-L)))+Lbar;