function [rar,yclean] = spm_rar (Z,p,m,verbose)
% Bayesian autoregressive modelling with zero-mean Gaussian mixture noise
% function [rar,yclean] = spm_rar (Z,p,m,verbose)
%
% Z          [N x 1] vector of data points
% p          Number of AR coefficients
% m          Number of mixture components (default=2)
% verbose    0/1 to printout inner workings (default=0)
%
% rar        Returned model 
% yclean     'Clean' data (ie. with outlier errors removed)
%
% -------------------------------------------------------
% The fields in rar are:
%
% p                The number of AR coefficients
% m                The number of components
% fm               The negative free energy
%
%                  In the field priors:
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
%
%                  In the field posts:
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% a_mean           AR parameters (posterior mean)
% a_cov            AR parameters (posterior cov)
% b_alpha,c_alpha  Gamma parameters for weight precisions
%
%                  Mean posterior values:
% pi               mixing coefficients (lambda/sum(lambda))
% variances        variances (1./(b.*c))
%
% gamma            the responsibilities of each noise component
%
% For details of algorithm see:
%
% S.J. Roberts and W.D. Penny. Variational Bayes for Generalised Autoregressive 
% models. IEEE Transactions on Signal Processing, 50(9):2245-2257, 2002
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_rar.m 1276 2008-03-28 18:29:19Z guillaume $

if nargin < 3 | isempty(m)
    m=2;
end

if nargin < 4 | isempty(verbose)
  verbose=0;
end

if (m==1)
  rar=spm_ar(Z,p);
  return;
end

N=length(Z);

% Initialise AR coefficients to maximum likelihood solution
Z=Z(:);
y=Z(p+1:end);
for i=1:p,
    x(:,i)=Z(p-i+1:end-i);
end
x=-x;
y2=y.^2;
xt=x';
N=size(x,1);
wun=ones(N,1);

a_mean = pinv(x)*y;
y_pred = x*a_mean;
err=y-y_pred;
v=mean((y-y_pred).^2);
a_cov = v*inv(x'*x);

% Set mixing priors
lambda_0=5;

% Setting to these values gives updates for mean_alpha
% v. close to evidence framework 
bishop_prior=1;
if bishop_prior
  b_alpha_prior=1000;
  c_alpha_prior=0.001;
  b_0=1000;
  c_0=0.001;
end
mean_alpha_prior=b_alpha_prior*c_alpha_prior;

% Cluster on absolute difference from mean
zmix=spm_kmeans1(abs(err-mean(err)),m);
% Posterior for mixers
lambda=100*zmix.pi;
zmean=[zmix.m].^2;
% Posterior for precisions
var_precision=1/std(zmean)^2;
for s=1:m,
      % Set so that b*c=precision and  b^2*c=var_precision
      precision=1/zmean(s);
      b(s)=var_precision/precision;
      c(s)=(precision^2)/var_precision;
end

% Initialise weight precision posterior
E_w=a_mean'*a_mean;
b_alpha=0.5*E_w+0.5*trace(a_cov)+(1/b_alpha_prior);
b_alpha=1/b_alpha;
c_alpha=0.5*p+c_alpha_prior;
mean_alpha=b_alpha*c_alpha;

if verbose
  disp('Init');
  disp('AR Coefficients');
  a_mean(:)'
  for s=1:m, 
    disp(sprintf('State %d mix=%1.2f var=%1.2f',s,lambda(s)/sum(lambda),1/(b(s)*c(s))));
  end
end

lik=[];
tol=0.0001;
max_loops=32;
WLOOPS=5;
for loops=1:max_loops,
    
    % E-step
    lambda_tot=sum(lambda);
    ypred=x*a_mean;
    ypred2=ypred.^2;
    y_err=sum(xt.*(a_cov*xt));
    tv=y2-2*ypred.*y+y_err'+ypred2;
    tv=tv';
    for s=1:m,
        log_tilde_pi(s)=psi(lambda(s))-psi(lambda_tot);
        log_tilde_beta(s)=psi(c(s))+log(b(s));
        tilde_pi(s)=exp(log_tilde_pi(s));
        tilde_beta(s)=exp(log_tilde_beta(s));
        mean_beta(s)=c(s)*b(s);
        tilde_var(s,:)=tv;
        gamma(s,:)=tilde_pi(s)*(tilde_beta(s)^0.5)*exp(-0.5*mean_beta(s)*tv);
    end
    gamma_n=sum(gamma);
    for s=1:m,
        if mean(gamma_n) > eps
            % If component still exists
            gamma(s,:)=gamma(s,:)./gamma_n;
        end
    end
    
    % M-step
    % Part I
    for s=1:m,
        pi_bar(s)=mean(gamma(s,:));
        N_bar(s)=N*pi_bar(s);
        mean_bar(s)=mean(gamma(s,:)'.*y);
        dg=diag(gamma(s,:));
        x_bar(s,:)=mean(dg*x);
        var_bar_mu(s)=mean(gamma(s,:).*tilde_var(s,:));
    end
    
    % CALCULATE THE FREE ENERGY
    avg_likelihood=sum(N_bar.*(log_tilde_pi+0.5*log_tilde_beta));
    fit=-0.5*N*sum(mean_beta.*var_bar_mu);
    % Ensure that 0 log 0 = 0
    ent_s=sum(sum(-gamma.*log(gamma+eps)));
    avg_likelihood=avg_likelihood+fit+ent_s;
    lambda_p=lambda_0*ones(1,m);
    kl_dir=spm_kl_dirichlet(lambda,lambda_p,log_tilde_pi);
    kl_gamm=0;
    for s=1:m,
        kl_gamm=kl_gamm+spm_kl_gamma(b(s),c(s),b_0,c_0);
    end
    kl_weights=spm_kl_normal(a_mean,a_cov,zeros(1,p),(1/mean_alpha)*eye(p));
    kl_alpha=spm_kl_gamma(b_alpha,c_alpha,b_alpha_prior,c_alpha_prior);
    fm= avg_likelihood - kl_dir - kl_gamm - kl_weights - kl_alpha;
    
    % Convergence criterion
    oldlik=lik;
    lik=fm;
    
    if (mod(loops-1,WLOOPS)==0)
        if (loops>1)
            if abs((lik-oldlik)/lik) < tol
                break;
            end
        end
    end
    
    % M-Step: Part II
    for s=1:m,
        % Mixers
        lambda(s)=N_bar(s)+lambda_0;
        % Precisions
        b(s)=1/(0.5*N*var_bar_mu(s)+1/b_0);
        c(s)=0.5*N_bar(s)+c_0;
        mean_beta(s)=c(s)*b(s);
    end
    
    if mod(loops,WLOOPS)==0
        % Weight precisions
        E_w=0.5*a_mean'*a_mean;
        b_alpha=E_w+0.5*trace(a_cov)+(1/b_alpha_prior);
        b_alpha=1/b_alpha;
        c_alpha=0.5*p+c_alpha_prior;
        mean_alpha=b_alpha*c_alpha;
        
        % AR coefficients
        cc=zeros(p,p);
        cw=zeros(p,1);
        for s=1:m,
            dg=diag(gamma(s,:));
            cc=cc+mean_beta(s)*x'*dg*x;
            cw=cw+mean_beta(s)*x'*dg*y;
        end
        cc=cc+mean_alpha*eye(p);
        a_cov=inv(cc);
        a_mean=a_cov*cw;
    end
    
    if verbose
        disp(sprintf('It=%d, L_AV =%1.2f, KL Mix=%1.2f, KL Prec=%1.2f, KL-AR=%1.2f, KL-alpha=%1.2f, Fm=%1.2f',loops,avg_likelihood,kl_dir,kl_gamm,kl_weights,kl_alpha,fm));
    end
    
end

% Put variables into data structure
rar.posts.a_mean=a_mean;
rar.posts.a_cov=a_cov;
rar.m=m;
rar.fm=fm;
rar.priors.lambda_0=lambda_0;
rar.priors.c_0=c_0;
rar.priors.b_0=b_0;

for k=1:m,
  rar.posts.lambda(k)=lambda(k);
end

for k=1:m,
  rar.posts.c(k)=c(k);
  rar.posts.b(k)=b(k);
end

for k=1:m,
  rar.pi(k)=lambda(k)/sum(lambda);
end
rar.variances=1./(b.*c);

% Pre-pad gamma with zeros to get original length time series
rar.gamma=[zeros(m,p),gamma];

% Get 'clean' data
if m > 1
    [tmp,outlier_class]=min(rar.pi);
    e=y-ypred;
    outlier_error=gamma(outlier_class,:)'.*e;
    yclean=ypred+e-outlier_error;
    yclean=[Z(1:p);yclean];  % Pre-padding
else
    yclean=Z;
end

