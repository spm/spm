function [x] = spm_inv_phi(y)
% inverse logistic function
% FORMAT [y] = spm_inv_phi(x)
%
% x   = log((y./(1 - y));
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_inv_phi.m 253 2005-10-13 15:31:34Z guillaume $

% apply
%--------------------------------------------------------------------------
x   = log(y./(1 - y));

