function [family,model] = spm_compare_families (lme,partition,names,ecp)
% Bayesian comparison of model families for group studies
% FORMAT [family,model] = spm_compare_families (lme,partition,names,ecp)
%
% lme           - array of log model evidences 
%                   rows: subjects
%                   columns: models (1..N)
% partition     - [1 x N] vector such that partition(m)=k signifies that
%                 model m belongs to family k (out of K) eg. [1 1 2 2 2 3 3]
% names         - cell array of K family names eg, {'fam1','fam2','fam3'}
% ecp           - compute exceedence probs ? (1 or 0, default=0)
%
% family        - family posterior  
%   .alpha0       prior counts
%   .alpha        posterior counts 
%   .r            expected values
%   .xp           exceedance probs
%
% model          - model posterior
%   .alpha0        prior counts
%   .alpha         posterior counts
%   .r             expected values
%
% This function assumes a uniform prior over model families (using a 
% prior count of unity for each family). It then 
% adjusts model priors accordingly, uses spm_BMS to get model posteriors,
% and computes family posteriors via aggregation.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_compare_families.m 3158 2009-05-28 16:28:27Z will $

if nargin < 4 | isempty(ecp)
    ecp=0;
end

% Number of samples for computing exceedance probs if K > 2
Nsamp=1e6;

% Number of models
N=length(partition);

% Number of families in partition
K=length(unique(partition));
family.alpha0=ones(1,K);

% Size of families 
for i=1:K,
    ind{i}=find(partition==i);
    fam_size(i)=length(ind{i});
end

% Set model priors to give uniform family prior
model.alpha0=zeros(1,N);
for i=1:K,
    model.alpha0(ind{i})=1/fam_size(i);
end

% Get model posterior
[alpha,exp_r,xp] = spm_BMS(lme, [], 0, 0, 0, model.alpha0);
model.alpha=alpha;
model.exp_r=exp_r;

% Get stats from family posterior
for i=1:K,
    family.alpha(i)=sum(model.alpha(ind{i}));
end
for i=1:K,
    family.exp_r(i)=family.alpha(i)/sum(family.alpha);
end
if ecp
    if N == 2
        % comparison of 2 families
        family.xp(1) = spm_Bcdf(0.5,family.alpha(2),family.alpha(1));
        family.xp(2) = spm_Bcdf(0.5,family.alpha(1),family.alpha(2));
    else
        % comparison of >2 families: use sampling approach
        family.xp = spm_dirichlet_exceedance(family.alpha,Nsamp);
    end
else
    family.xp=[];
end
