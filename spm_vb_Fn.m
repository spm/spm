function [F,L,KL] = spm_vb_Fn (Y,block)
% Compute voxel-wise contributions to model evidence
% FORMAT [F,L,KL] = spm_vb_Fn (Y,block)
%
% Y             [T x N] time series 
% block         data structure containing the following fields:
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_Fn.m 2451 2008-11-10 16:20:32Z lee $

T=block.T;
p=block.p;
k=block.k;
N=block.N;
X=block.X;

C2=block.C2;  

Bk = kron(diag(block.mean_alpha),block.Dw);
B=block.Hw*Bk*block.Hw';

if p>0
    if ~strcmp(block.priors.A,'Discrete')
        Jk = kron(diag(block.mean_beta),block.Da);
        J=block.Ha*Jk*block.Ha';
    end
end

% Get KL-alpha
KL.alpha=0;
for j = 1:k,
    KL.alpha=KL.alpha+spm_kl_gamma(block.b_alpha(j),block.c_alpha(j),block.b_alpha_prior(j),block.c_alpha_prior(j));
end

% Get KL-beta 
KL.beta=0;
if p > 0
    if strcmp(block.priors.A,'Discrete')
        for j=1:p,
            for s=1:block.priors.S,
                KL.beta=KL.beta+spm_kl_gamma(block.b_beta(j,s),block.c_beta(j,s),block.b_beta_prior(j,s),block.c_beta_prior(j,s));
            end
        end
    else
        for j = 1:p,
            KL.beta=KL.beta+spm_kl_gamma(block.b_beta(j),block.c_beta(j),block.b_beta_prior(j),block.c_beta_prior(j));
        end
    end
end

% Get average Likelihood, KL-Lambda and terms for KL-w, KL-a
for n=1:N,
    subblock_n           = [(n-1)*k+1:n*k];
    
    if p>0
        Gn=spm_vb_get_Gn(Y,block,n);
    else
        wc=block.w_cov{n};
        en=(Y(:,n)-X*block.w_mean(subblock_n,1));
        Gn=trace(wc*block.XTX)+en'*en;
    end
    
    L(n)=-0.5*block.mean_lambda(n)*Gn;
    L(n)=L(n) + 0.5*(T-p)*(psi(block.c_lambda(n)) + log(block.b_lambda(n)));
    L(n)=L(n)-0.5*block.C2/N;

    KL.lam(n)=spm_kl_gamma(block.b_lambda(n),block.c_lambda(n),block.b_lambda_prior(n),block.c_lambda_prior(n));
    
    KL_w1=-0.5*sum(log(block.mean_alpha))-0.5*log(det(block.w_cov{n}));
    KL_w1=KL_w1-0.5*k*block.log_det_Dw/N;
    KL_w2=0.5*trace(B(subblock_n,subblock_n)*block.w_cov{n});
    
    subblock_ni          = [1:N*k];
    subblock_ni(subblock_n) = [];
    Bnn               = B(subblock_n,subblock_n);
    Bni               = B(subblock_n,subblock_ni);
    Bin               = B(subblock_ni,subblock_n);
    
    w_mean=block.w_mean(subblock_n,1);
    KL_w3=0.5*w_mean'*Bnn*w_mean+0.5*w_mean'*Bni*block.w_mean(subblock_ni,1);
    KL_w3=KL_w3-0.5*k;
    KL.w(n)=KL_w1+KL_w2+KL_w3;
    
    if p>0
        asubblock_n  = [(n-1)*p+1:n*p];
        a_mean=block.a_mean(asubblock_n,1);
        if strcmp(block.priors.A,'Discrete')
            ibeta_n=diag(block.priors.gamma(n,:)*(1./block.mean_beta'));
            a_n=block.priors.gamma(n,:)*block.as';
            KL.a(n)=spm_kl_normal(a_mean,block.a_cov{n},a_n,ibeta_n);
%             a_mean
%             sqrt(block.a_cov{n})
%             a_n
%             sqrt(ibeta_n)
        else
            asubblock_ni = [1:N*p];
            asubblock_ni(asubblock_n)=[];
            Jnn = J(asubblock_n,asubblock_n);
            Jni = J(asubblock_n,asubblock_ni);
            KL_a1=-0.5*sum(log(block.mean_beta)) - 0.5*log(det(block.a_cov{n}));
            KL_a1=KL_a1-0.5*p*block.log_det_Da/N;
            KL_a2=0.5*trace(Jnn*block.a_cov{n});
            KL_a3=0.5*a_mean'*Jnn*a_mean+0.5*a_mean'*Jni*block.a_mean(asubblock_ni,1);
            KL_a3=KL_a3-0.5*p;
            KL.a(n)=KL_a1+KL_a2+KL_a3;
        end
    else
        KL.a(n)=0;
    end
end

F = L- KL.w  - KL.lam - KL.a - KL.alpha/N- KL.beta/N;
