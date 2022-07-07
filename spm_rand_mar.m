function [y] = spm_rand_mar(m,n,a)
% Generate random variates from an autoregressive process
% FORMAT [y] = spm_rand_mar(m,n,a)
% m   - time bins
% n   - variates
% a   - autoregression coefficients
%
% see also: spm_rand_power_law
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging
 

% create random process
%--------------------------------------------------------------------------
y  = spm_sqrtm(spm_Q(a,m))*randn(m,n);

