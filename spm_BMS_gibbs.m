function [exp_r,xp,r_samp] = spm_BMS_gibbs (lme, alpha0, Nsamp)
% Bayesian model selection for group studies using Gibbs sampling
% FORMAT [exp_r,xp,r_samp] = spm_BMS_gibbs (lme, alpha0, Nsamp)
%
% INPUT:
% lme      - array of log model evidences 
%              rows: subjects
%              columns: models (1..Nk)
% alpha0   - [1 x Nk] vector of prior model counts
% Nsamp    - number of samples (default: 1e6)
% 
% OUTPUT:
% exp_r   - [1 x  Nk] expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% r_samp  - [Nsamp x Nk] matrix of samples from posterior
% 
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_gibbs.m 3170 2009-06-01 12:03:31Z will $

if nargin < 3 || isempty(Nsamp)
    Nsamp = 1e6;
end

max_val = log(realmax('double'));
Ni      = size(lme,1);  % number of subjects
Nk      = size(lme,2);  % number of models

% prior observations
%--------------------------------------------------------------------------
if nargin < 3 || isempty(alpha0)
    alpha0  = ones(1,Nk);    
end
alpha0   = alpha0(:)';

% Initialise; sample r from prior
for k = 1:Nk
    r(:,k) = spm_gamrnd(alpha0(k),1);
end
sr= sum(r,2);
for k = 1:Nk
    r(:,k) = r(:,k)./sr;
end
    
% Gibbs sampling 
for samp=1:2*Nsamp,
    
    % Make this sparse later
    mod_vec=zeros(Ni,Nk);
    
    % Sample m's given y, r
    for i=1:Ni,
        % Pick a model for this subject
        u=exp(lme(i,:)+log(r))+eps;
        m=u/sum(u);
        mod=spm_multrnd(m,1);
        mod_vec(i,mod)=1;
    end
    
    % Sample r's given y, m
    beta=sum(mod_vec,1);
    alpha=alpha0+beta;
    for k = 1:Nk
        r(:,k) = spm_gamrnd(alpha(k),1);
    end
    sr= sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end

    % Only keep last Nsamp samples
    if samp > Nsamp
        r_samp(samp-Nsamp,:)=r;
    end
end

exp_r=mean(r_samp,1);
xp=[];
    