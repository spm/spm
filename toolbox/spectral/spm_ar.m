function [ar] = spm_ar (Z,p,verbose)
% Bayesian autoregressive modelling
% FORMAT [ar] = spm_ar (Z,p,verbose)
%
% y_pred (t) = -\sum_{i=1}^p a_i y (t-i) + e (t)
% Note the sign and ordering 
%
% The noise, e(t), is Gaussian
%
% Z             [N x 1] univariate time series 
% p             (scalar) order of model
% verbose       1=print out fitting progress (default=0)
%
% ar            data structure
% ----------------------------------
% ar.a_mean     AR coefficients
% ar.a_cov
% ar.mean_beta  error precision 
% ar.b_beta
% ar.c_beta
% ar.mean_alpha weight precision 
% ar.b_alpha
% ar.c_alpha
% ar.y          targets
% ar.y_pred     predictions
% ar.r2         proportion of variance explained
% ar.p          model order
% ar.fm         negative free energy
%
% For algorithmic details see:
%
% W.D. Penny and S.J. Roberts. Bayesian Methods for Autoregressive Models.
% In IEEE Workshop on Neural Networks for Signal Processing, Sydney Australia, 2000
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_ar.m 1276 2008-03-28 18:29:19Z guillaume $

if nargin < 2, 
   disp('spm_ar.m needs at least two arguments'); 
   return
end

if nargin < 3 | isempty(verbose)
    verbose=0;
end

Z=Z(:);
y=Z(p+1:end);
for i=1:p,
    x(:,i)=Z(p-i+1:end-i);
end

N=length(Z);
if abs(mean(Z)) > 3*(std(Z)/sqrt(N))
  % If mean is greater than 3SE away from 0
  % disp('Warning from vbar: mean subtracted from data');
  Z=Z-mean(Z);
end

% Initialise coefficients to maximum likelihood solution
% if p > 1
%     % In case columns of x are collinear
%     [ux,dx,vx]=svd(x);
%     ddx=diag(dx);
%     svd_tol=max(ddx)*eps*p;
%     rank_X=sum(ddx > svd_tol);
%     ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
%     ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
%     Xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';
%     X2=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)';
%     
%     a_mean= Xp*y;
%     y_pred= x*a_mean;
%     v=mean((y-y_pred).^2);
%     a_cov = v*X2;
%     xtx=X2;
% else
a_mean = pinv(x)*y;
y_pred = x*a_mean;
v=mean((y-y_pred).^2);
xtx=x'*x;
a_cov = v*inv(xtx);
    

% Setting to these values gives updates for mean_alpha, mean_beta
% approx to evidence framework 
b_alpha_prior=1000;
c_alpha_prior=0.001;
mean_alpha_prior=b_alpha_prior*c_alpha_prior;
b_beta_prior=1000;
c_beta_prior=0.001;

xty=x'*y;
xt=x';
yty=y'*y;

max_iters=32;
lik=[];
tol=0.0001;
for it=1:max_iters,

  E_w=a_mean'*a_mean;
  % Update weight precision
  b_alpha=0.5*E_w+0.5*trace(a_cov)+(1/b_alpha_prior);
  b_alpha=1/b_alpha;
  c_alpha=0.5*p+c_alpha_prior;
  mean_alpha=b_alpha*c_alpha;
  
  E_d_av=0.5*yty-a_mean'*xty;
  E_d_av=E_d_av+0.5*a_mean'*xtx*a_mean;
  E_d_av=E_d_av+0.5*trace(a_cov*xtx);

  % Update noise precision
  b_beta=1/(E_d_av+(1/b_beta_prior));
  c_beta=0.5*N+c_beta_prior;
  mean_beta=b_beta*c_beta;
  
  % Update weights
  a_cov=inv(mean_beta*xtx+mean_alpha*eye(p));
  a_mean=mean_beta*a_cov*xty;
  
  % Calculate f_m (negative free energy)
  l_av=0.5*N*(psi(c_beta)+log(b_beta))-0.5*N;
  kl_weights=spm_kl_normal(a_mean,a_cov,zeros(1,p),(1/mean_alpha)*eye(p));
  kl_alpha=spm_kl_gamma(b_alpha,c_alpha,b_alpha_prior,c_alpha_prior);
  kl_beta=spm_kl_gamma(b_beta,c_beta,b_beta_prior,c_beta_prior);
  f_m=l_av-kl_weights-kl_alpha-kl_beta;
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

% Now reverse sign of coefficients to keep format
% with ar_spec etc. (covariances will be unchanged)
ar.a_mean=-a_mean;

% Load up data structure
ar.y_pred = x*a_mean;
ar.y=y;
vy=std(y)^2;
ey=std(y-ar.y_pred)^2;
ar.r2=(vy-ey)/vy; 

ar.p=p;
ar.a_cov=a_cov;
ar.mean_beta=mean_beta;
ar.mean_alpha=mean_alpha;
ar.b_beta=b_beta;
ar.c_beta=c_beta;
ar.b_alpha=b_alpha;
ar.c_alpha=c_alpha;
ar.fm=f_m;
ar.l_av=l_av;
ar.iterations=it;
