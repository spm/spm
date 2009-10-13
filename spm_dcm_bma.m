function [theta] = spm_dcm_bma (post,subj,Nsamp,oddsr)
% Model-independent samples from DCM posterior  
% FORMAT [theta] = spm_dcm_bma (post,subj,Nsamp,oddsr)
%
% post      [Nd x M] vector of posterior model probabilities
%           If Nd>1 then inference is based on a RFX posterior p(r|Y)
% subj      subj(n).sess(s).model(m).fname: DCM filename
% Nsamp     Number of samples (default = 1e3)
% oddsr     posterior odds ratio for defining Occam's window (default=0, ie
%           all models used in average)
%
% theta     [Np x Nsamp] posterior density matrix. Parameter vector is of 
%           dimension Np and there are Nsamp samples
%
% This routine implements Bayesian averaging over models and subjects
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_dcm_bma.m 3461 2009-10-13 15:25:24Z maria $

if nargin < 3 | isempty(Nsamp)
    Nsamp=1e3;
end
if nargin < 4 | isempty(oddsr)
    oddsr=0;
end

Nsub=length(subj);
theta=[];

[Nd M]=size(post);
if Nd > 1 
    rfx=1;
    if Nsamp>Nd,
        error('Error in spm_dcm_bma: not enough samples');
    end
else
    rfx=0;
end

if rfx
    mean_post=mean(post);
    mp=max(mean_post);
    post_ind=find(mean_post>mp*oddsr);
    Nocc=length(post_ind);
    disp(' ');
    disp(sprintf('%d models in Occams window',Nocc));
    if Nocc==0, return; end
    for occ=1:Nocc,
        m=post_ind(occ);
        disp(sprintf('Model %d, <r|Y>=%1.2f',m,mean_post(m)));
    end
    
    % Renormalise post prob to Occam group
    post=post(:,post_ind);
    sp=sum(post,2);
    post=post./(sp*ones(1,Nocc));
else
    % Find models in Occam's window
    mp=max(post);
    post_ind=find(post>mp*oddsr);
    Nocc=length(post_ind);
    disp(' ');
    disp(sprintf('%d models in Occams window',Nocc));
    if Nocc==0, return; end
    for occ=1:Nocc,
        m=post_ind(occ);
        disp(sprintf('Model %d, p(m|Y)=%1.2f',m,post(m)));
    end
    
    % Renormalise post prob to Occam group
    post=post(post_ind);
    post=post/sum(post);
end

% Load DCM posteriors for models in Occam's window
for n=1:Nsub,
    for kk=1:Nocc,
        sel=post_ind(kk);
        load_str=subj(n).sess(1).model(sel).fname;
        load(load_str);
        subj(n).sess(1).model(kk).Ep=DCM.Ep;
        subj(n).sess(1).model(kk).Cp=full(DCM.Cp);
    end
end

% Pre-allocate sample arrays
Np=length(spm_vec(subj(n).sess(1).model(kk).Ep));
theta_all=zeros(Np,Nsub);

for i=1:Nsamp,
    % Pick a model
    if rfx
        m=spm_multrnd(post(i,:),1);
    else
        m=spm_multrnd(post,1);
    end
    
    % Pick parameters from model for each subject
    for n=1:Nsub,
        mu=subj(n).sess(1).model(m).Ep;
        mu=spm_vec(mu);
        sig=subj(n).sess(1).model(m).Cp;
        tmp=spm_samp_gauss (mu,sig,1)';
        theta_all(:,n)=tmp(:);
    end
    
    % Average over subjects
    if Nsub>1
        theta(:,i)=mean(theta_all,2);
    else
        theta(:,i)=theta_all;
    end
end

