function [x] = spm_inv_phi(y)
% inverse logistic function
% FORMAT [y] = spm_inv_phi(x)
%
% x   = log((y./(1 - y));
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_inv_phi.m 1143 2008-02-07 19:33:33Z spm $

% apply
%--------------------------------------------------------------------------
x   = log(y./(1 - y));

