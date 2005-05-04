function [F,L,KL] = spm_vb_Fn (Y,slice)
% Compute voxel-wise contributions to model evidence
% FORMAT [F,L,KL] = spm_vb_Fn (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure containing the following fields:
%
% F             [N x 1] vector where nth entry is unique contribution to 
%               model evidence from voxel n
% L             [N x 1] Average Likelihood
% KL.w          [N x 1] KL w - unique contribution
% KL.a          [N x 1] KL a - unique contribution
% KL.lam        [N x 1] KL Lambda
% KL.alpha      Scalar
% KL.beta       Scalar
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id: spm_vb_Fn.m 112 2005-05-04 18:20:52Z john $

T=slice.T;
p=slice.p;
k=slice.k;
N=slice.N;
X=slice.X;

C2=slice.C2;  

Bk = kron(diag(slice.mean_alpha),slice.Dw);
B=slice.Hw*Bk*slice.Hw';

if slice.p>0
    if ~strcmp(slice.priors.A,'Discrete')
        Jk = kron(diag(slice.mean_beta),slice.Da);
        J=slice.Ha*Jk*slice.Ha';
    end
end

% Get KL-alpha
KL.alpha=0;
for j = 1:k,
    KL.alpha=KL.alpha+spm_kl_gamma(slice.b_alpha(j),slice.c_alpha(j),slice.b_alpha_prior(j),slice.c_alpha_prior(j));
end

% Get KL-beta 
KL.beta=0;
if slice.p > 0
    if strcmp(slice.priors.A,'Discrete')
        for j=1:p,
            for s=1:slice.priors.S,
                KL.beta=KL.beta+spm_kl_gamma(slice.b_beta(j,s),slice.c_beta(j,s),slice.b_beta_prior(j,s),slice.c_beta_prior(j,s));
            end
        end
    else
        for j = 1:p,
            KL.beta=KL.beta+spm_kl_gamma(slice.b_beta(j),slice.c_beta(j),slice.b_beta_prior(j),slice.c_beta_prior(j));
        end
    end
end

% Get average Likelihood, KL-Lambda and terms for KL-w, KL-a
for n=1:slice.N,
    block_n           = [(n-1)*k+1:n*k];
    
    if slice.p>0
        Gn=spm_vb_get_Gn(Y,slice,n);
    else
        wc=slice.w_cov{n};
        en=(Y(:,n)-X*slice.w_mean(block_n,1));
        Gn=trace(wc*slice.XTX)+en'*en;
    end
    
    L(n)=-0.5*slice.mean_lambda(n)*Gn;
    L(n)=L(n) + 0.5*(T-p)*(spm_digamma(slice.c_lambda(n)) + log(slice.b_lambda(n)));
    L(n)=L(n)-0.5*slice.C2/slice.N;

    KL.lam(n)=spm_kl_gamma(slice.b_lambda(n),slice.c_lambda(n),slice.b_lambda_prior(n),slice.c_lambda_prior(n));
    
    KL_w1=-0.5*sum(log(slice.mean_alpha))-0.5*log(det(slice.w_cov{n}));
    KL_w1=KL_w1-0.5*k*slice.log_det_Dw/N;
    KL_w2=0.5*trace(B(block_n,block_n)*slice.w_cov{n});
    
    block_ni          = [1:N*k];
    block_ni(block_n) = [];
    Bnn               = B(block_n,block_n);
    Bni               = B(block_n,block_ni);
    Bin               = B(block_ni,block_n);
    
    w_mean=slice.w_mean(block_n,1);
    KL_w3=0.5*w_mean'*Bnn*w_mean+0.5*w_mean'*Bni*slice.w_mean(block_ni,1);
    KL_w3=KL_w3-0.5*k;
    KL.w(n)=KL_w1+KL_w2+KL_w3;
    
    if slice.p>0
        ablock_n  = [(n-1)*p+1:n*p];
        a_mean=slice.a_mean(ablock_n,1);
        if strcmp(slice.priors.A,'Discrete')
            ibeta_n=diag(slice.priors.gamma(n,:)*(1./slice.mean_beta'));
            a_n=slice.priors.gamma(n,:)*slice.as';
            KL.a(n)=spm_kl_normal(a_mean,slice.a_cov{n},a_n,ibeta_n);
%             a_mean
%             sqrt(slice.a_cov{n})
%             a_n
%             sqrt(ibeta_n)
%             keyboard
        else
            ablock_ni = [1:N*p];
            ablock_ni(ablock_n)=[];
            Jnn = J(ablock_n,ablock_n);
            Jni = J(ablock_n,ablock_ni);
            KL_a1=-0.5*sum(log(slice.mean_beta)) - 0.5*log(det(slice.a_cov{n}));
            KL_a1=KL_a1-0.5*p*slice.log_det_Da/N;
            KL_a2=0.5*trace(Jnn*slice.a_cov{n});
            KL_a3=0.5*a_mean'*Jnn*a_mean+0.5*a_mean'*Jni*slice.a_mean(ablock_ni,1);
            KL_a3=KL_a3-0.5*slice.p;
            KL.a(n)=KL_a1+KL_a2+KL_a3;
        end
    else
        KL.a(n)=0;
    end
end

F = L- KL.w  - KL.lam - KL.a - KL.alpha/slice.N- KL.beta/slice.N;
