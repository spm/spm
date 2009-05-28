function [family] = spm_compare_families (lme,partition,names,ecp)
% Bayesian comparison of model families for group studies
% FORMAT [family] = spm_compare_families (lme,partition,names,ecp)
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
%   .alpha        counts 
%   .r            expected values
%   .xp           exceedance probs
%
% This function assumes a uniform prior over model families. It then 
% adjusts model priors accordingly (to accomodate families containing
% unequal numbers of models), uses spm_BMS to get model posteriors,
% and then finally computes family posteriors via aggregation
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_compare_families.m 3155 2009-05-28 14:21:43Z will $

if nargin < 4 | isempty(ecp)
    ecp=0;
end

% Number of samples for computing exceedance probs if K > 2
Nsamp=1e6;

% Number of models
N=length(partition);

% Number of families in partition
K=length(unique(partition));

% Size of families 
for i=1:K,
    ind{i}=find(partition==i);
    fam_size(i)=length(ind{i});
end

% Get number of models in each family
fam_count=max(fam_size)./fam_size;

% Set model priors to give uniform family prior
alpha0=zeros(1,N);
for i=1:K,
    alpha0(ind{i})=fam_count(i);
end

% Get model posterior
[alpha,exp_r,xp] = spm_BMS(lme, [], 0, 0, 0, alpha0);

% Get stats from family posterior
for i=1:K,
    family.alpha(i)=sum(alpha(ind{i}));
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
