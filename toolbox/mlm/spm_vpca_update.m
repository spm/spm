function [pca,c] = spm_vpca_update (T,S,pca,c,m)
% Update VPCA parameters
% FORMAT [pca,c] = spm_vpca_update (T,S,pca,c,m)
%
% T     [d x N] matrix containing N d-dimensional data vectors
%       The nth data sample, t_n, is nth column of T
% S
% pca   data structure (see eg. spm_vpca.m)
% c     information about single component
% m     cluster number (used for mixtures of VPCA model)
%
% pca,c updated info
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vpca_update.m 5962 2014-04-17 12:47:43Z spm $

q=c(m).q;
d=pca.d;
N=pca.N;
if isfield(pca,'update_mean')
    update_mean=pca.update_mean;
else
    update_mean=1;
end

% Update latent variables
c(m).sxn2=zeros(q,q);
c(m).Sigma_x=inv(eye(q)+pca.mean_tau*c(m).avg_WtW);
for n=1:N,
    c(m).M_x(:,n)=pca.mean_tau*c(m).Sigma_x*c(m).M_w'*(T(:,n)-c(m).mean_mu);
    c(m).xn2(:,:,n)=c(m).Sigma_x+c(m).M_x(:,n)*c(m).M_x(:,n)';
    c(m).sxn2=c(m).sxn2+S(m,n)*c(m).xn2(:,:,n);
end

% Update factors 
c(m).Sigma_w = inv(diag(c(m).mean_alpha)+pca.mean_tau*c(m).sxn2);
for k=1:d,
    obs_err=T(k,:)-c(m).mean_mu(k)*ones(1,N);
    obs_err=S(m,:).*obs_err;
    x_err=c(m).M_x*obs_err';
    mwk=pca.mean_tau*c(m).Sigma_w*x_err;
    c(m).M_w(k,:) = mwk'; 
end
c(m).avg_WtW=c(m).M_w'*c(m).M_w+d*c(m).Sigma_w;


% Update mean
if update_mean
    c(m).Sigma_mu=1/(pca.beta+sum(S(m,:))*pca.mean_tau);
    t_err=zeros(d,1);
    for n=1:N,
        t_err=t_err+S(m,n)*(T(:,n)-c(m).M_w*c(m).M_x(:,n));
    end
    c(m).mean_mu=pca.mean_tau*c(m).Sigma_mu*t_err;
end

% Update alphas - shrinkage priors on factor columns
pca.qa_alpha=pca.a_alpha+d/2;
for i=1:q,
    c(m).qb_alpha(i)=pca.b_alpha+0.5*c(m).avg_WtW(i,i);
    c(m).mean_alpha(i)=pca.qa_alpha/c(m).qb_alpha(i);
end

c(m).sum_quad_exp=0;
T3=0.5*spm_logdet(c(m).Sigma_x);
quad_term=c(m).mean_mu'*c(m).mean_mu+trace(c(m).Sigma_mu);
trSx=trace(c(m).Sigma_x);
mutW=2*c(m).mean_mu'*c(m).M_w;
for n=1:N,
    term1=T(:,n)'*T(:,n)+quad_term;
    term2=trace(c(m).avg_WtW*c(m).xn2(:,:,n));
    term3=mutW*c(m).M_x(:,n);
    term4=-2*T(:,n)'*c(m).M_w*c(m).M_x(:,n)-2*T(:,n)'*c(m).mean_mu;
    c(m).sum_quad_exp=c(m).sum_quad_exp+S(m,n)*(term1+term2+term3+term4);
    c(m).T1b(n)=-0.5*(c(m).M_x(:,n)'*c(m).M_x(:,n)+trSx);
    T2=-0.5*pca.mean_tau*(term1+term2+term3+term4);
    pca.log_QS(m,n)=c(m).T1b(n)+T2+T3;
end

if pca.M==1
    % Update tau - precision of observation noise
    pca.qa_tau=pca.a_tau+0.5*N*d;
    pca.qb_tau=pca.b_tau+0.5*c(m).sum_quad_exp;
    pca.mean_tau=pca.qa_tau/pca.qb_tau;
end
