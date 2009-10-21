function [theta, Nocc] = spm_dcm_bma (post,post_indx,subj,Nsamp,oddsr)
% Model-independent samples from DCM posterior  
% FORMAT [theta, Nocc] = spm_dcm_bma (post,post_indx,subj,Nsamp,oddsr)
%
% post      [Ni x M] vector of posterior model probabilities
%           If Ni>1 then inference is based on subject-specific RFX posterior 
% post_indx models to use in BMA (position of models in subj sctructure)
% subj      subj(n).sess(s).model(m).fname: DCM filename
% Nsamp     Number of samples (default = 1e3)
% oddsr     posterior odds ratio for defining Occam's window (default=0, ie
%           all models used in average)
%
% theta     [Np x Nsamp] posterior density matrix. Parameter vector is of 
%           dimension Np and there are Nsamp samples
% Nocc      Number of models in Occam's window
%
%           For RFX BMA, different subject can have different models in
%           Occam's window (and different numbers of models in Occam's
%           window)
%
% This routine implements Bayesian averaging over models and subjects
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_dcm_bma.m 3493 2009-10-21 14:55:53Z will $

if nargin < 3 | isempty(Nsamp)
    Nsamp=1e3;
end
if nargin < 4 | isempty(oddsr)
    oddsr=0;
end

Nsub=length(subj);
theta=[];

[Ni M]=size(post);
if Ni > 1 
    rfx=1;
else
    rfx=0;
end

if rfx
    for i=1:Ni,
        mp=max(post(i,:));
        post_ind{i}=find(post(i,:)>mp*oddsr);
        Nocc(i)=length(post_ind{i});
        disp(' ');
        disp(sprintf('Subject %d has %d models in Occams window',i,Nocc(i)));
        if Nocc(i)==0,
            return;
        end
        
        for occ=1:Nocc(i),
            m=post_ind{i}(occ);
            disp(sprintf('Model %d, <p(m|Y>=%1.2f',m,post(i,m)));
        end
        
        % Renormalise post prob to Occam group
        renorm(i).post=post(i,post_ind{i});
        sp=sum(renorm(i).post,2);
        renorm(i).post=renorm(i).post./(sp*ones(1,Nocc(i)));
        
        % Load DCM posteriors for models in Occam's window
        for kk=1:Nocc(i),
            sel=post_indx(post_ind{i}(kk));
            load_str=subj(i).sess(1).model(sel).fname;
            load(load_str);
            params(i).model(kk).Ep=DCM.Ep;
            params(i).model(kk).Cp=full(DCM.Cp);
        end
    end
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
    
    % Load DCM posteriors for models in Occam's window
    for n=1:Nsub,
        for kk=1:Nocc,
            sel=post_indx(post_ind(kk));
            load_str=subj(n).sess(1).model(sel).fname;
            load(load_str);
            params(n).model(kk).Ep=DCM.Ep;
            params(n).model(kk).Cp=full(DCM.Cp);
        end
    end
end

% Pre-allocate sample arrays
Np=length(spm_vec(params(1).model(1).Ep));
theta_all=zeros(Np,Nsub);

for i=1:Nsamp,
    % Pick a model
    if ~rfx
        m=spm_multrnd(post,1);
    end
    % Pick parameters from model for each subject
    for n=1:Nsub,
        if rfx
            m=spm_multrnd(renorm(n).post,1);
        end
        mu=params(n).model(m).Ep;
        mu=spm_vec(mu);
        sig=params(n).model(m).Cp;
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

