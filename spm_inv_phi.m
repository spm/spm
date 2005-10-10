function [x] = spm_inv_phi(y)
% inverse logistic function
% FORMAT [y] = spm_inv_phi(x)
%
% x   = log((y./(1 - y));
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% apply
%--------------------------------------------------------------------------
x   = log(y./(1 - y));

