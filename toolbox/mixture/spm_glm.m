function [rglm] = spm_glm (y,X,alpha,verbose)
% Fit a Bayesian GLM 
% FORMAT [rglm] = spm_glm (y,X,alpha,verbose)
%
% This function is called by spm_robust_glm if m==1
% 
% y          [N x 1] data vector
% X          [N x p] design matrix
% alpha      [1 x 1] weight precision (default=0.001)
% verbose    0/1 to printout inner workings (default=0)
%
% rglm       Returned model 
%
% -------------------------------------------------------
% The fields in rglm are:
%
% m                The number of error components
% fm               The negative free energy
%
%                  In the field priors:
%
% b_0,c_0          Gamma parameters for precisions
%
%                  In the field posts:
%
% b,c              Gamma parameters for precisions
% w_mean           Mean estimated regression coefficients
% w_cov            Covariance of regression coefficients
%
%                  Mean posterior values:
% variances        variances (1./(b.*c))
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_glm.m 1276 2008-03-28 18:29:19Z guillaume $

if nargin < 4 | isempty(verbose)
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

% Prior for beta
b_beta_prior=1000;
c_beta_prior=0.001;

% Initialise weight precision
if nargin < 3 | isempty(alpha)
    mean_alpha=0.001;
else
    mean_alpha=alpha;
end

xtx=X'*X;
xty=X'*y;
xt=X';
yty=y'*y;

max_iters=100;
lik=[];
tol=0.0001;
for it=1:max_iters,
  
  E_d_av=0.5*yty-w_mean'*xty;
  E_d_av=E_d_av+0.5*w_mean'*xtx*w_mean;
  E_d_av=E_d_av+0.5*trace(w_cov*xtx);

  % Update noise precision
  b_beta=1/(E_d_av+(1/b_beta_prior));
  c_beta=0.5*N+c_beta_prior;
  mean_beta=b_beta*c_beta;
  
  % Update weights
  w_cov=inv(mean_beta*xtx+mean_alpha*eye(p));
  w_mean=mean_beta*w_cov*xty;
  
  % Calculate f_m (negative free energy)
  l_av=0.5*N*(psi(c_beta)+log(b_beta))-mean_beta*E_d_av;
  kl_weights=spm_kl_normald(w_mean,w_cov,zeros(1,p),(1/mean_alpha)*eye(p));
  kl_beta=spm_kl_gamma(b_beta,c_beta,b_beta_prior,c_beta_prior);
  f_m=l_av-kl_weights-kl_beta;
  
  if verbose
      disp(sprintf('Iteration %d L_av=%1.3f KL_w=%1.3f KL_alpha=%1.3f KL_beta=%1.3f Fm=%1.3f',it,l_av,kl_weights,kl_alpha,kl_beta,f_m));
  end

  % Convergence criterion
  oldlik=lik;
  lik=f_m;
  
  if (it<=2)
    likbase=lik;
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase))
    break;
  end;

end

% Put variables into data structure
rglm.its=it;
rglm.m=1;
rglm.fm=f_m;
rglm.kl_gamm=kl_beta;
rglm.kl_w=kl_weights;
rglm.priors.c_0=c_beta_prior;
rglm.priors.b_0=b_beta_prior;
rglm.w=w_mean;
rglm.posts.w_mean=w_mean;
rglm.posts.w_cov=w_cov;
rglm.posts.c=c_beta;
rglm.posts.b=b_beta;
rglm.variances=1./(b_beta.*c_beta);

