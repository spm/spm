function [logbf] = spm_mlm_posthoc (mlm,c,a)
% Post-hoc model comparison of multivariate linear models
% FORMAT [logbf] = spm_mlm_posthoc (mlm,c,a)
%
% mlm          MLM data structure - see spm_mlm_bayes.m
%              This contains eg. the [p x d] posterior mean regression  
%              coefficient matrix mlm.wmean. 
%
% c            [k x p*d] contrast matrix defining k-dimensional subspace
% a            hypothesized value (zeros(k,1) by default)
%
% The contrast matrix and hypothesized value define the reduced model. 
% The contrast is applied to the vectorised parameters w = vec(mlm.wmean)
%             
% The Bayes Factor in favour of the alternative hypothesis over the null
% is computed using a Savage-Dickey ratio (the probability of the
% hypothesized value under the prior versus its probability under the
% posterior)
%
% bf = p(c*w=a|mlm)/p(c*w=a|Y,mlm)              
%
% logbf        Log Bayes Factor
%___________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mlm_posthoc.m 4651 2012-02-09 16:03:39Z will $

wpost=mlm.wmean(:);

post_mean=c*wpost;
k=size(post_mean,1);
prior_mean=zeros(k,1);

post_cov=c*mlm.wcov*c';

% Get prior covariance matrix
prec_prior=size(mlm.wcov,1);
for g=1:mlm.prior.groups,
    prec_prior=prec_prior+mlm.prior.group(g).mean_alpha*mlm.prior.group(g).index;
end
w_prior_cov=diag(1./prec_prior);
prior_cov=c*w_prior_cov*c';

if nargin < 3 | isempty (a)
    a=zeros(k,1);
end

% Enforce symmetry - to avoid numerical issues
prior_cov=(prior_cov+prior_cov')/2;
post_cov=(post_cov+post_cov')/2;

% For now just use spm_mvNpdf. This can be made more efficient !
prior_prob=spm_mvNpdf(a,prior_mean,prior_cov);
post_prob=spm_mvNpdf(a,post_mean,post_cov);

logbf=log(prior_prob)-log(post_prob);