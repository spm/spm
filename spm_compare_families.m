function [family,model] = spm_compare_families (lme,partition,names,Nsamp,prior)
% Bayesian comparison of model families for group studies using Gibbs sampling
% FORMAT [family,model] = spm_compare_families (lme,partition,names,Nsamp,prior)
%
% lme           - array of log model evidences 
%                   rows: subjects
%                   columns: models (1..N)
% partition     - [1 x N] vector such that partition(m)=k signifies that
%                 model m belongs to family k (out of K) eg. [1 1 2 2 2 3 3]
% names         - cell array of K family names eg, {'fam1','fam2','fam3'}
% Nsamp         - Number of samples to get 
% prior         - 'F-unity' alpha0=1 for each family (default)
%               - 'M-unity' alpha0=1 for each model (not advised)
%
% family        - family posterior  
%   .alpha0       prior counts
%   .exp_r        expected value of r
%   .rsamp        samples from posterior
%   .xp           exceedance probs
%
% model          - model posterior
%   .alpha0        prior counts
%   .exp_r         expected value of r
%
% This function computes model/family posteriors using samples from spm_BMS_gibbs.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_compare_families.m 3236 2009-06-29 17:08:03Z will $

if nargin < 4 | isempty(Nsamp)
    Nsamp=1e4;
end

if nargin < 5 | isempty(prior)
    prior='F-unity';
end

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

% Set model priors 
switch prior,
    case 'F-unity',
        for i=1:K,
            model.alpha0(ind{i})=1/fam_size(i);
        end
    case 'M-unity',
        model.alpha0=ones(1,N);
    otherwise
        disp('Error in spm_compare_families:Unknown prior');
end

% Get model posterior
[exp_r,tmp,r_samp]=spm_BMS_gibbs(lme,model.alpha0,Nsamp);
model.exp_r=exp_r;

% Get stats from family posterior
for i=1:K,
    ri=r_samp(:,ind{i});
    family.r_samp(:,i)=sum(ri,2);
    family.exp_r(i)=mean(family.r_samp(:,i));
end

% Exceedence probs
xp = zeros(1,K);
r=family.r_samp;
[y,j]=max(r,[],2);
tmp=histc(j,1:K)';
family.xp=tmp/Nsamp;
