function [mix] = spm_mix (y,m,verbose)
% Fit a multivariate Gaussian Mixture model using VB
% FORMAT [mix] = spm_mix (y,m,verbose)
%
% y          [N x d] data matrix containing N samples of d-dim data
% m          Number of mixture components
% verbose    Set to 1 to see evolution of free energy, 0 otherwise
%            (default=1)
%
% mix        Returned model
%
%--------------------------------------------------------------------------
% The fields in mix are:
%
% m                The number of components
% fm               The negative free energy. This decomposes into
%                  fm=acc-kl_proportions-kl_covs-kl_centres
%
% acc              model accuracy
% kl_proportions   complexity penalty for cluster proportions
% kl_covs          complexity penalty for cluster covariances
% kl_centres       complexity penalty for cluster centres
%
% Fields:
%
% lambda           Post mixers, q(pi|D) = D(lambda)
% gamma            [m x N] matrix of belonging probabilities
% state(s).a       Post precisions, q(Gamma|D)=W(a,B)
% state(s).B   
% state(s).C       Post covariance
% state(s).m       Post mean, q(mu|D)=N(m_s,beta_s Gamma_s)
% state(s).beta
% state(s).prior   Estimated mixing proportions
%
%                  In the field prior:
%
% lambda_0         Prior mixers, p(pi) = D(lambda_0)
% a_0,B_0          Prior precisions, p(Gamma)=W(a_0,B_0)
% m_0,beta_0       Prior means, p(mu)=N(m_0,beta_0 Gamma_s)
%
%__________________________________________________________________________
% Copyright (C) 2007-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mix.m 5962 2014-04-17 12:47:43Z spm $

% This code implements the algorithm in:
%
% See Attias, H. (2000) A Variational
% Bayesian framework for Graphical Models, NIPS.
%
% See also:
%
% W.D. Penny (2001) Variational Bayes for d-dimensional Gaussian mixture models 
% Wellcome Department of Imaging Neuroscience, University College London. 
%
% for Negative Free Energy expression and validation.


if nargin < 3
    verbose=1;
end
[N,d]=size(y);

% Put model info into data structure format
mix.nin=d;
mix.m=m;

% PRIORS
% Set mixing priors
lambda_0=1;
% Set mean priors
m_0=mean(y);
beta_0=1;
% Set precision priors
% These settings set the prior covariance matrix to be
% equal to the 0.01*d*identity matrix - if original data is
% zmuv, this is reasonable
a_0=d; 
B_0=0.01*d*eye(d);
%B_0=d*eye(d);

% Also, store in data structure ready for output
mix.prior.lambda_0=lambda_0;
mix.prior.m_0=m_0;
mix.prior.beta_0=beta_0;
mix.prior.a_0=a_0;
mix.prior.B_0=B_0;

sd=0;
for i=1:d,
  sd=sd+psi((a_0+1-i)/2);
end
log_tilde_gamma_0=sd-log(det(B_0))+d*log(2);
   
if (m==1)
    % For a single component we just have a Gaussian model
    
    % Set posterios
    state(1).a=N+a_0;
    state(1).beta=N+beta_0;
    state(1).m=m_0;
    Cy=cov(y);
    state(1).B=N*Cy+B_0;
    
    fm=spm_lg_gamma (d,0.5*state(1).a);
    fm=fm-spm_lg_gamma (d,0.5*a_0);
    fm=fm-0.5*d*N*log(pi);
    fm=fm+0.5*d*(log(beta_0)-log(state(1).beta));
    fm=fm+0.5*a_0*log(det(B_0)); 
    fm=fm-0.5*state(1).a*log(det(state(1).B)); 
    
    mix.m=m;
    mix.fm=fm;
    state(1).C=state(1).B/state(1).a;
    mix.state=state;
    
    mix.priors=1;
    return;
end

% Run kmeans on the data
% to initialise posterior 
[priors,means,covs] = spm_kmeans(y,m);

% Add pseudo-counts to ML priors
priors=priors+(1/N);
lambda=priors;
for s=1:m,
  state(s).m=means(s,:);
  state(s).beta=priors(s)*N+beta_0;
  state(s).a=priors(s)*N+a_0;
  Cov=covs(:,:,s);
  % If component cov is rank deficient set to prior
  if rank(Cov) < size(Cov,1)
    state(s).B=B_0;
  else
    state(s).B=priors(s)*N*Cov;
  end
  N_bar(s)=priors(s)*N;
end

% Start algorithm
lik=[];
tol=0.0001;
max_loops=32;
for loops=1:max_loops,
    
    % E-step
    lambda_tot=sum(lambda);
    for s=1:m,
        state(s).bar_gamma=state(s).a*inv(state(s).B);
        log_tilde_pi(s)=psi(lambda(s))-psi(lambda_tot);
        sd=0;
        for i=1:d,
            sd=sd+psi((state(s).a+1-i)/2);
        end
        log_tilde_gamma(s)=sd-log(det(state(s).B))+d*log(2);
        
        tilde_pi(s)=exp(log_tilde_pi(s));
        tilde_gamma(s)=exp(log_tilde_gamma(s));
        for n=1:N,
            gamma(s,n)=tilde_pi(s)*tilde_gamma(s)^0.5;
            dy=(y(n,:)-state(s).m);
            gamma(s,n)=gamma(s,n)*(exp(-0.5*dy*state(s).bar_gamma*dy')+eps)*exp(-d/(2*state(s).beta));
        end
    end
    
    gamma_n=sum(gamma);
    for s=1:m,
        if mean(gamma_n) > eps
            % If component still exists
            gamma(s,:)=gamma(s,:)./gamma_n;
        end
    end
    
    % M-step
    
    % Part-I
    for s=1:m,
        pi_bar(s)=mean(gamma(s,:))+eps;
        N_bar(s)=N*pi_bar(s)+eps;
        bar_mu(s,:)=(1/N_bar(s))*sum(gamma(s,:)'*ones(1,d).*y);
        
    end
    
    % get weighted means and covariances
    for s=1:m,
        state(s).bar_sigma=zeros(d,d);
        for n=1:N,
            dy=y(n,:)-bar_mu(s,:);
            state(s).bar_sigma=state(s).bar_sigma+gamma(s,n).*(dy'*dy);
        end
        state(s).bar_sigma=(state(s).bar_sigma)/N_bar(s);
        
        if isnan(state(s).bar_sigma)
            state(s).bar_sigma
        end
    end
    
    % Now compute the free energy 
    f1=-spm_kl_dirichlet(lambda,lambda_0*ones(1,m),log_tilde_pi);
    for s=1:m,
        f2(s)=-spm_kl_wishart(state(s).a,state(s).B,a_0,B_0);
        
        % KL-method for computing f3(s)
        Cs=state(s).B/(state(s).beta*state(s).a);
        C0=state(s).B/(beta_0*state(s).a);
        f3(s)=-spm_kl_normal(state(s).m,Cs,m_0,C0);
        
        % Check f3(s)
        check_f3=-0.5*d*(log(state(s).beta)-log(beta_0)+(beta_0/state(s).beta)-1);
        dm=state(s).m-m_0;
        check_f3=check_f3-0.5*dm*beta_0*state(s).a*inv(state(s).B)*dm';
        f3(s)=check_f3;
        
        f4(s)=N_bar(s)*log_tilde_pi(s)-sum(gamma(s,:).*log(gamma(s,:)+eps));
        LaB=0;
        for i=1:d,
            LaB=LaB+psi((state(s).a+1-i)/2);
        end  
        LaB=LaB+d*log(2)-log(det(state(s).B));
        f5(s)=0.5*N_bar(s)*(-d*log(2*pi)+LaB-trace(state(s).bar_gamma*state(s).bar_sigma)-(d/state(s).beta));
        fkl(s)=f2(s)+f3(s)+f4(s)+f5(s);
        fkl_adj(s)=f2(s)+f3(s)+(f4(s)+f5(s))/N_bar(s);
    end
    fm=f1+sum(f2)+sum(f3)+sum(f4)+sum(f5);
 
    acc=sum(f4)+sum(f5);
    kl_proportions=-f1;
    kl_covs=-sum(f2);
    kl_centres=-sum(f3);
    
    if verbose
        disp(sprintf('Iter=%d, F1=%1.2f, F2=%1.2f, F3=%1.2f, F4=%1.2f, F5=%1.2f, Fm=%1.2f',loops,f1,sum(f2),sum(f3),sum(f4),sum(f5),fm));
    end 
 
    % Convergence criterion
    oldlik=lik;
    lik=fm;
    if (loops>1)
        if abs((lik-oldlik)/lik) < tol
            break;
        end
    end;
    
    % Part-II
    
    for s=1:m,
        % Posterior mixing coefficients
        lambda(s)=N_bar(s)+lambda_0;
        % Posterior means
        state(s).m=(N_bar(s)*bar_mu(s,:)+beta_0*m_0)/(N_bar(s)+beta_0);
        state(s).beta=N_bar(s)+beta_0;
        % Posterior precisions
        state(s).a=N_bar(s)+a_0;
        dy=bar_mu(s,:)-m_0;
        state(s).B=N_bar(s)*state(s).bar_sigma + (N_bar(s)*beta_0*dy'*dy)/(N_bar(s)+beta_0)+B_0;
    end
end

% Put variables into data structure
mix.m=m;
mix.fm=fm;
mix.acc=acc;
mix.kl_proportions=kl_proportions;
mix.kl_covs=kl_covs;
mix.kl_centres=kl_centres;
mix.state=state;
mix.lambda=lambda;
mix.gamma=gamma;

% Put info into data structure
for j=1:m,
    mix.state(j).prior=pi_bar(j);
    mix.state(j).C=mix.state(j).B/mix.state(j).a;
end
