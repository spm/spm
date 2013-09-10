function [y] = spm_rand_mar(m,n,a)
% generates random variates from an autoregressive process
% FORMAT [y] = spm_rand_mar(m,n,a)
% m   - time bins
% n   - variates
% a   - autoregression coefficients
%
% see also: spm_rand_power_law
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_rand_mar.m 5633 2013-09-10 13:58:03Z karl $
 

% create random process
%--------------------------------------------------------------------------
y  = spm_sqrtm(spm_Q(a,m))*randn(m,n);

