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
% $Id: spm_BMS_gibbs.m 3252 2009-07-06 18:07:53Z guillaume $

if nargin < 3 || isempty(Nsamp)
    Nsamp = 1e3;
end

Ni      = size(lme,1);  % number of subjects
Nk      = size(lme,2);  % number of models

% prior observations
%--------------------------------------------------------------------------
if nargin < 3 || isempty(alpha0)
    alpha0  = ones(1,Nk);    
end
alpha0   = alpha0(:)';

% Initialise; sample r from prior
r  = zeros(1,Nk);
for k = 1:Nk
    r(:,k) = spm_gamrnd(alpha0(k),1);
end
sr = sum(r,2);
for k = 1:Nk
    r(:,k) = r(:,k)./sr;
end

% Subtract subject means
lme=lme-mean(lme,2)*ones(1,Nk);

% Gibbs sampling 
r_samp = zeros(Nsamp,Nk);
for samp=1:2*Nsamp
    
    mod_vec=sparse(Ni,Nk);
    % Sample m's given y, r
    for i=1:Ni
        % Pick a model for this subject
        u=exp(lme(i,:)+log(r))+eps;
        m=u/sum(u);
        modnum=spm_multrnd(m,1);
        mod_vec(i,modnum)=1;
    end
    
    % Sample r's given y, m
    beta=sum(mod_vec,1);
    alpha=alpha0+beta;
    for k = 1:Nk
        r(:,k) = spm_gamrnd(alpha(k),1);
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end

    % Only keep last Nsamp samples
    if samp > Nsamp
        r_samp(samp-Nsamp,:)=r;
    end
    
    if mod(samp,1e4)==0
        disp(sprintf('%d samples out of %d',samp,2*Nsamp));
    end
    
end

exp_r = mean(r_samp,1);
xp    = [];
    