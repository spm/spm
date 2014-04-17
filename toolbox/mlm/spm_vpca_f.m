function [Fm] = spm_vpca_f (pca,c) 
% Compute free energy of VPCA model
% FORMAT [Fm] = spm_vpca_f (pca,c)
%
% pca   data structure (see eg. spm_vpca.m)
% c     information about single component
%
% Fm    negative free energy of model
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vpca_f.m 5962 2014-04-17 12:47:43Z spm $


d=pca.d;
q=c(1).q;
pca.q=q;

N=pca.N;
pca.sum_quad_exp=c(1).sum_quad_exp;
pca.T1b=c(1).T1b;
pca.Sigma_x=c(1).Sigma_x;
pca.qb_alpha=c(1).qb_alpha;
pca.mean_alpha=c(1).mean_alpha;
pca.M_w=c(1).M_w;
pca.Sigma_w=c(1).Sigma_w;
pca.mean_mu=c(1).mean_mu;
pca.Sigma_mu=c(1).Sigma_mu;

l_av1 = psi(pca.qa_tau) - log(pca.qb_tau); % Term 1
l_av1 = 0.5*pca.N*pca.d*l_av1 - 0.5*pca.N*pca.d*log(2*pi);
l_av1 = l_av1 - 0.5 * pca.mean_tau * pca.sum_quad_exp;

I_T=l_av1;

l_av2 = - 0.5*pca.N*pca.q * log (2*pi) + sum(pca.T1b);

l_av3=0.5*pca.N*pca.q*(1+log(2*pi))+0.5*pca.N*spm_logdet(pca.Sigma_x);

I_X=l_av2+l_av3;

l_av=l_av1+l_av2+l_av3;

W_term=0;
W_term=-0.5*q*d*log(2*pi);
dig_qa=psi(pca.qa_alpha);
for i=1:q,
    exp_log_alpha_i= dig_qa - log(pca.qb_alpha(i));
    W_term=W_term+0.5*d*exp_log_alpha_i;
end
MAlpha=diag(pca.mean_alpha);
for k=1:d,
    W_term=W_term-0.5*pca.M_w(k,:)*MAlpha*pca.M_w(k,:)';
end
W_term=W_term-0.5*d*trace(MAlpha*pca.Sigma_w);

% Entropy of q(W)
for k=1:d,
    W_term=W_term+0.5*q*(1+log(2*pi))+0.5*spm_logdet(pca.Sigma_w);
end
I_W=W_term;

kl_alpha=0;
for i=1:q,
    kl_alpha=kl_alpha+spm_kl_gamma (1/pca.qb_alpha(i), pca.qa_alpha, 1/pca.b_alpha, pca.a_alpha);
end
    
kl_mu=0;
for k=1:d,
    kl_mu = kl_mu+spm_kl_normald (pca.mean_mu(k),pca.Sigma_mu,0,1/pca.beta);
end

kl_tau = spm_kl_gamma (1/pca.qb_tau, pca.qa_tau, 1/pca.b_tau, pca.a_tau);

I_alpha=-kl_alpha;
I_mu=-kl_mu;
I_tau=-kl_tau;

kl_sum = kl_alpha + kl_mu + kl_tau;
Fm = l_av + W_term - kl_sum;
