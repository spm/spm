function [y] = spm_phi(x)
% Logistic function
% FORMAT [y] = spm_phi(x)
%
% y   = 1./(1 + exp(-x))
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% apply
%--------------------------------------------------------------------------
y   = 1./(1 + exp(-x));
