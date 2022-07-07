function [y] = spm_phi_dot(x)
% Return the derivative of the logistic function
% FORMAT [y] = spm_phi_dot(x)
% see spm_phi and spm_inv_phi
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% apply
%--------------------------------------------------------------------------
u   = exp(-x);
y   = 1./(1+u).^2.*u;
