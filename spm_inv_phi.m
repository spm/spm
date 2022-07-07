function [x] = spm_inv_phi(y)
% Inverse logistic function
% FORMAT [y] = spm_inv_phi(x)
%
% x   = log((y./(1 - y))
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% apply
%--------------------------------------------------------------------------
x = log(y./(1 - y));
