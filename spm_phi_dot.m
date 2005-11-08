function [y] = spm_phi_dot(x)
% returns the derivative of the logistic function
% FORMAT [y] = spm_phi_dot(x)
% see spm_phi and spm_inv_phi
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% apply
%--------------------------------------------------------------------------
u   = exp(-x);
y   = 1./(1+u).^2.*u;
