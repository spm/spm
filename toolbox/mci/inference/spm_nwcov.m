function [M] = spm_nwcov (M)
% Get second moments of Normal-Wishart
% FORMAT [M] = spm_nwcov (M)
%
% .mean_prior_cov    Prior covariance of mean
% .sample_prior_cov  Prior covariance of samples
% .mean_post_cov     Posterior covariance of mean
% .sample_pred_cov   Predictive covariance of samples
%
% The latter quantity is also the covariance of the predictive density
% The marginal distributions of the mean and of the samples 
% are multivariate-T, not Gaussian.
%
% See J. Bernardo and A. Smith (2000) 
% Bayesian Theory, Wiley (page 435)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_nwcov.m 6548 2015-09-11 12:39:47Z will $

% Get prior covariances
M.mean_prior_cov=M.B0/(M.n0*(M.a0-1));

alpha=2*M.a0-M.P+1;
w_s=(1+1/M.n0)/(0.5*alpha-1);
M.sample_prior_cov=w_s*M.B0;

% Get posterior covariances
M.mean_post_cov=M.BN/(M.nN*(M.aN-1));

alpha=2*M.aN-M.P+1;
w_s=(1+1/M.nN)/(0.5*alpha-1);
M.sample_pred_cov=w_s*M.BN;