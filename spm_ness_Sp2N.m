function [m,C] = spm_ness_Sp2N(Sp,n,K)
% Convert polynomial potential parameters into a Gaussian density
% FORMAT [m,C] = spm_ness_Sp2N(Sp,[n,K])
%--------------------------------------------------------------------------
% Sp - Polynomial coefficients or parameters of log density
% n  - Dimensionality of state space
% K  - Order of polynomial expansion (K = 3 corresponds to quadratic)
%
% m  - (Gaussian) mean
% C  - (Gaussian) covariance
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 3, K = 3; end
if nargin < 2, n = 3; end

% get sufficient statistics of Gaussian density
%--------------------------------------------------------------------------
[m,C] = spm_ness_cond(n,K,Sp);
