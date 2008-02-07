function [y] = spm_phi_dot(x)
% returns the derivative of the logistic function
% FORMAT [y] = spm_phi_dot(x)
% see spm_phi and spm_inv_phi
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_phi_dot.m 1143 2008-02-07 19:33:33Z spm $

% apply
%--------------------------------------------------------------------------
u   = exp(-x);
y   = 1./(1+u).^2.*u;
