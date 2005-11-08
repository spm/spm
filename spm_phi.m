function [y] = spm_phi(x)
% logistic function
% FORMAT [y] = spm_phi(x)
%
% y   = 1./(1 + exp(-x));
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% apply
%---------------------------------------------------------------------------
y   = 1./(1 + exp(-x));

