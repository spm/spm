function [y] = spm_phi(x)
% logistic function
% FORMAT [y] = spm_phi(x)
%
% y   = 1./(1 + exp(-x));
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_phi.m 289 2005-11-10 17:15:04Z guillaume $

% apply
%---------------------------------------------------------------------------
y   = 1./(1 + exp(-x));

