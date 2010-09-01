function spm_modify_resms(VResMS)
% Modification of artefactually low residual mean squared error estimates
% FORMAT spm_modify_resms(VResMS) % Updates VResMS in place.
% 
% This function implements a heuristic shrinkage estimate, inspired by the
% formal procedures of James & Stein [1] and Ledoit & Wolf [2].
% The default method of operation (defaults.stats.resms.method = 'max') is
% motivated by the prior belief that three orders of magnitude difference
% (tunable with defaults.stats.resms.amount) in variance should allow 
% reasonable variation, while avoiding problematically low variances.
%
% Refences:
% [1] James, W. & Stein, C. (1961) Estimation with quadratic loss.
%     Fourth Berkeley Symposium on Mathematical Statistics and Probability.
%     http://tinyurl.com/james1961estimation
% [2] Ledoit, O. & Wolf, M. (2004) A well-conditioned estimator for
%     large-dimensional covariance matrices. 
%     Journal of Multivariate Analysis, 88:365-411
%     http://dx.doi.org/10.1016/S0047-259X(03)00096-4
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Ged Ridgway
% $Id: spm_modify_resms.m 4058 2010-09-01 14:26:34Z ged $

% defaults.stats.resms.method allows for future variations (e.g. use of
% cortical mask to avoid influence of variance within deep MEG sources).

method = spm_get_defaults('stats.resms.method');
amount = spm_get_defaults('stats.resms.amount');

ResMS = spm_read_vols(VResMS);

switch method
    case 'off'
        return
    case 'max'
        ResMS = ResMS + amount * max(ResMS(isfinite(ResMS)));
    otherwise
        error('Unrecognised shrinkage method %s', method)
end

spm_write_vol(VResMS, ResMS);
