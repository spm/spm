function [F,Lav, KL] = spm_vb_F (Y,slice)
% Compute lower bound on evidence, F, for VB GLMAR models
% FORMAT [F,Lav, KL] = spm_vb_F (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure containing the following fields:
%
% F             Lower bound on model evidence, F
% Lav           Average Likelihood
% KL            Kullback-Liebler Divergences with fields
%               .w, .alpha, .beta, .Lambda, .a
%
% This function implements equation 18 in paper VB4
%
% Will Penny & Nelson Trujillo-Barreto $Id$

if slice.verbose
    disp('Updating F');
end

T=slice.T;
p=slice.p;
k=slice.k;
N=slice.N;
X=slice.X;

C2=slice.C2;  

Bk = kron(diag(slice.mean_alpha),slice.D);
B=slice.Hw*Bk*slice.Hw';

if slice.p>0
    Jk = kron(diag(slice.mean_beta),slice.D);
    J=slice.Ha*Jk*slice.Ha';
end

tr_B_qcov=0;
log_det_qcov=0;

if slice.p>0
    tr_J_acov=0;
    log_det_acov=0;
    KL_a=0;
end

% Get average Likelihood, KL-Lambda and terms for KL-w
KL_lambda=0;
C1=0;
Lav_term=0;
for n=1:slice.N,
    block_n           = [(n-1)*k+1:n*k];
    ablock_n          = [(n-1)*p+1:n*p];
    if slice.p > 0
        G(n,1)=spm_vb_get_Gn(Y,slice,n);
        tr_J_acov=tr_J_acov+trace(J(ablock_n,ablock_n)*slice.a_cov{n});
        log_det_acov=log_det_acov+log(det(slice.a_cov{n}));
    else
        wc=slice.w_cov{n};
        en=(Y(:,n)-X*slice.w_mean(block_n,1));
        Gn=trace(wc*slice.XTX)+en'*en;
        Lav_term=Lav_term+slice.mean_lambda(n)*Gn;
    end

    C1  = C1 + spm_digamma(slice.c_lambda(n)) + log(slice.b_lambda(n));
    KL_lambda=KL_lambda+spm_kl_gamma(slice.b_lambda(n),slice.c_lambda(n),slice.b_lambda_prior(n),slice.c_lambda_prior(n));
    
    tr_B_qcov=tr_B_qcov+trace(B(block_n,block_n)*slice.w_cov{n});
    log_det_qcov=log_det_qcov+log(det(slice.w_cov{n}));
end
if slice.p > 0
    Lav_term=slice.mean_lambda.'*G;
end
Lav = ((T-p)*C1 - Lav_term - C2)./2;

% Get KL-alpha
KL_alpha=0;
log_det_alphas=0;
for j = 1:k,
    KL_alpha=KL_alpha+spm_kl_gamma(slice.b_alpha(j),slice.c_alpha(j),slice.b_alpha_prior(j),slice.c_alpha_prior(j));
    log_det_alphas=log_det_alphas+log(slice.mean_alpha(j));
end
term1=-0.5*N*log_det_alphas;

% Get KL-beta 
if slice.p > 0
    KL_beta=0;
    log_det_betas=0;
    for j = 1:p,
        KL_beta=KL_beta+spm_kl_gamma(slice.b_beta(j),slice.c_beta(j),slice.b_beta_prior(j),slice.c_beta_prior(j));
        log_det_betas=log_det_betas+log(slice.mean_beta(j));
    end
    beta_term1=-0.5*N*log_det_betas;
end

% Get KL-w
try 
    % We only log determinant of spatial precision matrix
    % if we are comparing different spatial priors
    slice.log_det_D;
catch
    [vvv,ddd]=eig(full(slice.D));
    slice.log_det_D=sum(log(diag(ddd)));
end
term1=term1-0.5*k*slice.log_det_D;
KL_w=term1-0.5*log_det_qcov+0.5*tr_B_qcov+0.5*slice.w_mean'*B*slice.w_mean-0.5*N*k;

% Get KL-a and add up terms to get F
if slice.p > 0
    beta_term1=beta_term1-0.5*p*slice.log_det_D;
    KL_a=beta_term1-0.5*log_det_acov+0.5*tr_J_acov+0.5*slice.a_mean'*J*slice.a_mean-0.5*N*p;
    F = Lav - (KL_w+KL_alpha+KL_lambda+KL_a+KL_beta);
else
    F = Lav - (KL_w+KL_alpha+KL_lambda);
    KL_a=0;KL_beta=0;
end

KL.w=KL_w;
KL.alpha=KL_alpha;
KL.beta=KL_beta;
KL.Lambda=KL_lambda;
KL.a=KL_a;
