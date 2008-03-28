function [rglm,yclean] = spm_rglm (y,X,m,priors,verbose)
% Fit a Robust GLM 
% FORMAT [rglm,yclean] = spm_rglm (y,X,m,priors,verbose)
%
% The noise is modelled with a Mixture of Zero-Mean Gaussians 
%
% y          [N x 1] data vector
% X          [N x p] design matrix
% m          Number of mixture components
% priors     .alpha      [1 x 1] weight precision (default=0.001)
%            .mean_err   [m x 1] vector of mean error SD
%            .std_err    [m x 1] vector of dev of error SD
% verbose    0/1 to printout inner workings (default=0)
%
% rglm       Returned model 
% yclean     'Clean' data
%
% -------------------------------------------------------
% The fields in rglm are:
%
% m                The number of error components
% fm               The negative free energy
% loops            Number of iterations used
%
%                  In the field priors:
%
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
%
%                  In the field posts:
%
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% w_mean           Mean estimated regression coefficients
% w_cov            Covariance of regression coefficients
% pi               mixing coefficients (lambda/sum(lambda))
% variances        variances (1./(b.*c))
% gamma            the responsilities of each noise component
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_rglm.m 1276 2008-03-28 18:29:19Z guillaume $

if nargin < 4 | isempty(priors)
    mean_alpha=0.001;
    % Set variance priors
    b_0=1000*ones(m,1);
    c_0=0.001*ones(m,1);
else
    mean_alpha=priors.alpha;
    for mm=1:m,
        b_0(mm)=(priors.std_err(mm)^2)/priors.mean_err(mm);
        c_0(mm)=priors.mean_err(mm)/b_0(mm);
    end
end

if (m==1)
    rglm=spm_glm(y,X,mean_alpha);
    return
end

if nargin < 5 | isempty(verbose)
    verbose=0;
end

y=y(:);
N=length(y);
p=size(X,2);
spX=issparse(X);

% Initialise regression coefficients to maximum likelihood solution
if spX
    iX=inv(X'*X);
    w_mean=iX*X'*y;
else
    iX=inv(X'*X);
    w_mean = pinv(X)*y;
end
y_pred = full(X*w_mean);
err=y-y_pred;
v=mean((y-y_pred).^2);
w_cov = v*iX;

% Set mixing priors
lambda_0=5;

% Cluster on absolute difference from mean
zmix=spm_kmeans1(abs(err-mean(err)),m);
% Posterior for mixers
lambda=100*zmix.pi+lambda_0;
zmean=[zmix.m].^2;
% Posterior for precisions
var_precision=1/std(zmean)^2;
for s=1:m,
      % Set so that b*c=precision and  b^2*c=var_precision
      precision=1/zmean(s);
      b(s)=var_precision/precision;
      c(s)=(precision^2)/var_precision;
end

if verbose
    disp('Init');
    for s=1:m, 
        disp(sprintf('State %d mix=%1.2f var=%1.2f',s,lambda(s)/sum(lambda),1/(b(s)*c(s))));
    end
end

lik=[];
tol=0.0001;
max_loops=100;
WLOOPS=5;
for loops=1:max_loops,
    
    % E-step
    lambda_tot=sum(lambda);
    ypred=full(X*w_mean);
    ypred2=ypred.^2;
    y2=y.^2;
    xt=X';
    y_err=full(sum(xt.*(w_cov*xt)));
    
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
        dg=sparse(1:N,1:N,gamma(s,:));
        x_bar(s,:)=mean(dg*X);
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
        kl_gamm=kl_gamm+spm_kl_gamma(b(s),c(s),b_0(s),c_0(s));
    end
    kl_weights=spm_kl_normald(w_mean,w_cov,zeros(1,p),(1/mean_alpha)*eye(p));
    fm= avg_likelihood - kl_dir - kl_gamm - kl_weights;
    
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
        b(s)=1/(0.5*N*var_bar_mu(s)+1/b_0(s));
        c(s)=0.5*N_bar(s)+c_0(s);
        mean_beta(s)=c(s)*b(s);
    end
    
    if mod(loops,WLOOPS)==0
        % Regression coefficients
        cc=sparse(p,p);
        cw=zeros(p,1);
        for s=1:m,
            dg=sparse(1:N,1:N,gamma(s,:));
            cc=cc+mean_beta(s)*X'*dg*X;
            cw=cw+mean_beta(s)*X'*dg*y;
        end
        cc=cc+mean_alpha*speye(p);
        w_cov=inv(cc);
        w_mean=w_cov*cw;
    end
    
    if verbose
        disp(sprintf('It=%d, L_AV =%1.2f, KL Mix=%1.2f, KL Prec=%1.2f, KL-W=%1.2f, Fm=%1.2f',loops,avg_likelihood,kl_dir,kl_gamm,kl_weights,fm));
        disp(' ');
        disp(sprintf('Iteration number=%d',loops));
        for s=1:m,
            var=1/(b(s)*c(s));
            disp(sprintf('State %d pi_bar=%1.2f var=%1.2f',s,pi_bar(s),var));
        end
    end
    
    
end

% Put variables into data structure
rglm.m=m;
rglm.fm=fm;
rglm.kl_dir=kl_dir;
rglm.kl_gamm=kl_gamm;
rglm.kl_w=kl_weights;
rglm.mean_alpha=mean_alpha;
rglm.priors.lambda_0=lambda_0;
rglm.priors.c_0=c_0;
rglm.priors.b_0=b_0;
for k=1:m,
  rglm.posts.lambda(k)=lambda(k);
end
for k=1:m,
  rglm.posts.c(k)=c(k);
  rglm.posts.b(k)=b(k);
end
rglm.posts.w_mean=w_mean;
rglm.posts.w_cov=w_cov;

for k=1:m,
  rglm.posts.pi(k)=lambda(k)/sum(lambda);
end
rglm.posts.variances=1./(b.*c);
rglm.posts.gamma=gamma;
rglm.loops=loops;

% Get 'clean' data
if m > 1
    [tmp,outlier_class]=min(rglm.posts.pi);
    e=y-ypred;
    outlier_error=rglm.posts.gamma(outlier_class,:)'.*e;
    yclean=ypred+e-outlier_error;
else
    yclean=y;
end
