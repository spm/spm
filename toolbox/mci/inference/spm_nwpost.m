function [M] = spm_nwpost (M,w)
% Get posterior distribution over m,Lambda
% FORMAT [M] = spm_nwpost (M,w)
%
% M     M.prior - params of Normal-Wishart prior
% w     Multivariate data samples
%
% M     M.post - params of Normal-Wishart posterior
%
% Bernardo and Smith, Bayesian Theory, 2000 (p.441)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_nwpost.m 6548 2015-09-11 12:39:47Z will $

mw=mean(w,2);
Sw=cov(w',1);

N=size(w,2);

prior=M.prior;
P=prior.P;
a0=prior.a;
B0=prior.B;
beta0=prior.beta;
m0=prior.m;

post.beta=beta0+N;
post.m=(beta0*m0+N*mw)/post.beta;

post.a=a0+N/2;
post.B=B0+0.5*N*Sw+0.5*(beta0*N/post.beta)*(mw-m0)*(mw-m0)';
post.a=a0+N/2;

% Quantities for predictive density (over new samples)
post.mu_w=post.m;
w_s=(post.beta/(post.beta+1))*(post.a-0.5*(P-1));
post.Lambda_w=w_s*inv(post.B);
post.v_w=2*post.a-P+1;

% Wrap up
post.P=P;
post.N=N;
M.post=post;