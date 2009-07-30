function [family,model] = spm_compare_families (lme,family)
% Bayesian comparison of model families for group studies 
% FORMAT [family,model] = spm_compare_families (lme,family)
%
% INPUT:
%
% lme           - array of log model evidences 
%                   rows: subjects
%                   columns: models (1..N)
%
% family        - data structure containing family definition and inference parameters:
%                  .infer='RFX' or 'FFX' (default)
%                  .partition  [1 x N] vector such that partition(m)=k signifies that
%                              model m belongs to family k (out of K) eg. [1 1 2 2 2 3 3]
%                  .names      cell array of K family names eg, {'fam1','fam2','fam3'}
%                  .Nsamp      RFX only: Number of samples to get (default=1e4)
%                  .prior      RFX only: 'F-unity' alpha0=1 for each family (default)
%                              or 'M-unity' alpha0=1 for each model (not advised)
%
% OUTPUT:
%
% family        - RFX only:  
%                   .alpha0       prior counts 
%                   .exp_r        expected value of r
%                   .rsamp        samples from posterior
%                   .xp           exceedance probs
%                - FFX only: 
%                   .prior        family priors
%                   .post         family posteriors
%
% model          - RFX only: 
%                   .alpha0        prior counts
%                   .exp_r         expected value of r
%                - FFX only: 
%                   .subj_lme      log model ev without subject effects
%                   .prior         model priors
%                   .like          model likelihoods
%                   .posts         model posteriors
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_compare_families.m 3298 2009-07-30 15:30:29Z will $

try
    infer=family.infer;
catch
    disp('Error in spm_compare_families: inference method not specified');
    return
end

try
    partition=family.partition;
catch
    disp('Error in spm_compare_families: partition not specified');
    return
end

try
    names=family.names;
catch
    disp('Error in spm_compare_families: names not specified');
    return
end

if strcmp(infer,'RFX')
    try
        Nsamp=family.Nsamp;
    catch
        Nsamp=1e4;
        family.Nsamp=Nsamp;
    end
    
    try
        prior=family.prior;
    catch
        prior='F-unity';
        family.prior='F-unity';
    end
end

% Number of models
N=length(partition);

% Number of families in partition
K=length(unique(partition));

% Size of families 
for i=1:K,
    ind{i}=find(partition==i);
    fam_size(i)=length(ind{i});
end

if strcmp(infer,'FFX')
    
    % Family priors
    for i=1:K,
        family.prior(i)=1/K;
    end
    
    % Model priors
    for i=1:N,
        model.prior(i)=1/fam_size(partition(i));
    end
    
    % Model likelihoods
    lme=lme-mean(lme,2)*ones(1,N); % Subtract subject effects
    model.subj_lme=lme;
    model.like=sum(lme,1);
    model.like=exp(model.like);
    
    % Model posterior
    num=model.prior.*model.like;
    model.post=num/sum(num);
    
    % Family posterior
    for i=1:K,
        family.post(i)=sum(model.post(ind{i}));
    end
    
    return;
end
    
% Set model priors 
switch prior,
    case 'F-unity',
        for i=1:K,
            model.alpha0(ind{i})=1/fam_size(i);
        end
        family.alpha0=ones(1,K);
    case 'M-unity',
        model.alpha0=ones(1,N);
        for i=1:K,
            family.alpha0(i)=fam_size(i);
        end
    otherwise
        disp('Error in spm_compare_families:Unknown prior');
end

% Get model posterior
[exp_r,xp,r_samp]=spm_BMS_gibbs(lme,model.alpha0,Nsamp);
model.exp_r=exp_r;
model.xp=xp;

% Get stats from family posterior
for i=1:K,
    ri=r_samp(:,ind{i});
    family.r_samp(:,i)=sum(ri,2);
    family.exp_r(i)=mean(family.r_samp(:,i));
end

% Family exceedence probs
xp = zeros(1,K);
r=family.r_samp;
[y,j]=max(r,[],2);
tmp=histc(j,1:K)';
family.xp=tmp/Nsamp;
