function [j] = spm_mci_gprior_deriv (Pr,M)
% Gradient of Log Gaussian prior
% FORMAT [j] = spm_mci_gprior_deriv (Pr,M)
%
% Pr        parameters (vectorised and in M.V subspace)
% M         model structure
%
% j         gradient of log Gaussian prior
%__________________________________________________________________________

% Will Penny and Biswa Sengupta
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

Pr=Pr(:);

% Parameters errors in subspace
e = Pr;

% Gradient
j = - e'*M.ipC;
