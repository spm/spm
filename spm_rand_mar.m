function [y] = spm_rand_mar(m,n,a)
% generates random variates from an autoregressive process
% FORMAT [y] = spm_rand_mar(m,n,a)
% m   - time bins
% n   - variates
% a   - autoregression coefficients
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_phase_shuffle.m 3334 2009-08-25 16:13:38Z karl $
 

% create random process
%--------------------------------------------------------------------------
y  = spm_sqrtm(spm_Q(a,m))*randn(m,n);

