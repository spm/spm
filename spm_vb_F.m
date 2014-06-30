function [F,Lav,KL] = spm_vb_F(Y,block)
% Compute lower bound on evidence, F, for VB-GLM-AR models
% FORMAT [F,Lav,KL] = spm_vb_F(Y,block)
%
% Y             [T x N] time series 
% block         data structure (see spm_vb_glmar)
%
% F             Lower bound on model evidence, F
% Lav           Average Likelihood
% KL            Kullback-Liebler Divergences with fields
%               .w, .alpha, .beta, .Lambda, .a
%
% This function implements equation 18 in paper VB4.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_F.m 6079 2014-06-30 18:25:37Z spm $

if block.verbose
    disp('Updating F');
end

T = block.T;
p = block.p;
k = block.k;
N = block.N;
X = block.X;

C2 = block.C2;  

Bk = kron(diag(block.mean_alpha),block.Dw);
B  = block.Hw*Bk*block.Hw';

if p>0
    if ~strcmp(block.priors.A,'Discrete')
        Jk = kron(diag(block.mean_beta),block.Da);
        J  = block.Ha*Jk*block.Ha';
    end
end

tr_B_qcov     = 0;
log_det_qcov  = 0;

if p>0
    tr_J_acov = 0;
    log_det_acov = 0;
    KL_a      = 0;
end

% Get average Likelihood, KL-Lambda and terms for KL-w
KL_lambda = 0;
C1 = 0;
Lav_term = 0;
for n=1:N
    subblock_n  = [(n-1)*k+1:n*k];
    asubblock_n = [(n-1)*p+1:n*p];
    if p > 0
        G(n,1)  = spm_vb_get_Gn(Y,block,n);
        if ~strcmp(block.priors.A,'Discrete')
            tr_J_acov    = tr_J_acov+trace(J(asubblock_n,asubblock_n)*block.a_cov{n});
            log_det_acov = log_det_acov+log(det(block.a_cov{n}));
        end
    else
        wc       = block.w_cov{n};
        en       = (Y(:,n)-X*block.w_mean(subblock_n,1));
        Gn       = trace(wc*block.XTX)+en'*en;
        Lav_term = Lav_term+block.mean_lambda(n)*Gn;
    end

    C1  = C1 + psi(block.c_lambda(n)) + log(block.b_lambda(n));
    KL_lambda = KL_lambda + spm_kl_gamma(block.b_lambda(n),block.c_lambda(n),block.b_lambda_prior(n),block.c_lambda_prior(n));
    
    tr_B_qcov = tr_B_qcov+trace(B(subblock_n,subblock_n)*block.w_cov{n});
    log_det_qcov = log_det_qcov+log(det(block.w_cov{n}));
end
if p > 0
    Lav_term=block.mean_lambda.'*G;
end
Lav = ((T-p)*C1 - Lav_term - C2)./2;

%-Get KL-alpha
%--------------------------------------------------------------------------
KL_alpha = 0;
log_det_alphas = 0;
for j=1:k
    KL_alpha = KL_alpha + spm_kl_gamma(block.b_alpha(j),block.c_alpha(j),block.b_alpha_prior(j),block.c_alpha_prior(j));
    log_det_alphas = log_det_alphas+log(block.mean_alpha(j));
end
term1 = -0.5*N*log_det_alphas;

%-Get KL-beta
%--------------------------------------------------------------------------
if p > 0
    KL_beta = 0;
    if strcmp(block.priors.A,'Discrete')
        for j=1:p
            for s=1:block.priors.S
                KL_beta = KL_beta + spm_kl_gamma(block.b_beta(j,s),block.c_beta(j,s),block.b_beta_prior(j,s),block.c_beta_prior(j,s));
            end
        end
    else
        log_det_betas = 0;
        for j=1:p
            KL_beta = KL_beta + spm_kl_gamma(block.b_beta(j),block.c_beta(j),block.b_beta_prior(j),block.c_beta_prior(j));
            log_det_betas = log_det_betas + log(block.mean_beta(j));
        end
        beta_term1 = -0.5*N*log_det_betas;
    end
end

% Get KL-w
term1 = term1 -0.5*k*block.log_det_Dw;
KL_w  = term1 -0.5*log_det_qcov+0.5*tr_B_qcov+0.5*block.w_mean'*B*block.w_mean-0.5*N*k;

% Get KL-a and add up terms to get F
if p > 0
    if strcmp(block.priors.A,'Discrete')
        KL_a = 0;
        for n=1:N
            subblock_n = [(n-1)*p+1:n*p]; % Index for AR coeffs
            a_mean     = block.a_mean(subblock_n,1);
            ibeta_n    = diag(block.priors.gamma(n,:)*(1./block.mean_beta'));
            a_n        = block.priors.gamma(n,:)*block.as';
            KL_a       = KL_a + spm_kl_normal(a_mean,block.a_cov{n},a_n,ibeta_n);
        end
    else 
        beta_term1     = beta_term1 -0.5*p*block.log_det_Da;
        KL_a           = beta_term1 -0.5*log_det_acov+0.5*tr_J_acov+0.5*block.a_mean'*J*block.a_mean-0.5*N*p;
    end
    F       = Lav - (KL_w + KL_alpha + KL_lambda + KL_a + KL_beta);
else
    F       = Lav - (KL_w + KL_alpha + KL_lambda);
    KL_a    = 0;
    KL_beta = 0;
end

KL.w      = KL_w;
KL.alpha  = KL_alpha;
KL.beta   = KL_beta;
KL.Lambda = KL_lambda;
KL.a      = KL_a;
