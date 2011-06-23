function [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp)
% Bayesian model selection for group studies using Gibbs sampling
% FORMAT [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp)
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
% g_post  - [Ni x Nk] matrix of posterior probabilities with 
%           g_post(i,k) being post prob that subj i used model k
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_gibbs.m 4376 2011-06-23 15:03:16Z will $

if nargin < 3 || isempty(Nsamp)
    Nsamp = 1e4;
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

% Ensure all log model evidence differences are now within machine range
max_val = log(realmax('double'));
for i=1:Ni,
    for k = 1:Nk,
        lme(i,k) = sign(lme(i,k)) * min(max_val,abs(lme(i,k)));
    end
end
        
% Gibbs sampling 
r_samp = zeros(Nsamp,Nk);
g_post = zeros(Ni,Nk);

for samp=1:2*Nsamp
    
    mod_vec=sparse(Ni,Nk);
    % Sample m's given y, r
    for i=1:Ni
        % Pick a model for this subject
        u=exp(lme(i,:)+log(r))+eps;
        g=u/sum(u);
        gmat(i,:)=g;
        modnum=spm_multrnd(g,1);
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
        g_post=g_post+gmat;
    end
    
    if mod(samp,1e4)==0
        disp(sprintf('%d samples out of %d',samp,2*Nsamp));
    end
    
end
g_post=g_post/Nsamp;

% Posterior mean
exp_r = mean(r_samp,1);

% Exceedence probs
xp = zeros(1,Nk);
[y,j]=max(r_samp,[],2);
tmp=histc(j,1:Nk)';
xp=tmp/Nsamp;


    