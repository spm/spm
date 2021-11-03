function [m,C] = spm_ness_Sp2N(Sp,n,K)
% converts polynomial potential parameters into a Gaussian density
% FORMAT [m,C] = spm_ness_Sp2N(Sp,[n,K])
%--------------------------------------------------------------------------
% Sp - Polynomial coefficients or parameters of log density
% n  - Dimensionality of state space
% K  - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% m  - (Gaussian) mean
% C  - (Gaussian) covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_cond.m 8097 2021-04-24 20:28:27Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3, K = 3; end
if nargin < 2, n = 3; end

% get sufficient statistics of Gaussian density
%--------------------------------------------------------------------------
[m,C] = spm_ness_cond(n,K,Sp);

return

