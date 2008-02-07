function [y] = spm_phi(x)
% logistic function
% FORMAT [y] = spm_phi(x)
%
% y   = 1./(1 + exp(-x));
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_phi.m 1143 2008-02-07 19:33:33Z spm $

% apply
%---------------------------------------------------------------------------
y   = 1./(1 + exp(-x));

