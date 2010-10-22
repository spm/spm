function [mar,y,y_pred] = spm_mar (X,p,prior,verbose)
% Bayesian Multivariate Autoregressive Modelling
% FORMAT [mar,y,y_pred] = spm_mar (X,p,prior,verbose)
%
% Matrix of AR coefficients are in form
% x_t = -a_1 x_t-1 - a_2 x_t-2 + ...... - a_p x_t-p
% where a_k is a d-by-d matrix of coefficients at lag k and x_t-k's are 
% vectors of a d-variate time series.
%
% X              T-by-d matrix containing d-variate time series
% p              Order of MAR model
% prior          Prior on MAR coefficients (see marprior.m)
% verbose        1 to print out iteration details, 0 otherwise (default=0)
%
% mar.lag(k).a   AR coefficient matrix at lag k
% mar.noise_cov  Estimated noise covariance
% mar.fm         Free energy of model
% mar.wmean      MAR coefficients stored in a matrix
% y              Target values
% y_pred         Predicted values
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar.m 4096 2010-10-22 19:40:34Z karl $

if nargin < 4 || isempty(verbose)
    verbose=0;
end

d=size(X,2);    % dimension of time series
N=size(X,1);    % length of time series

% Embedding of multiple time series 
% giving x=[(x1(t-p) x2(t-p) .. xd(t-p)) (x1(t-p+1) x2(t-p+1)..xd(t-p+1)) ...
%           (x1(t-1) x2(t-1) .. xd(t-1))] on each row

x=[];
for i=1:p,
  tmpx=X(i:N-p+i,:);
  x=[x,tmpx];
end

% Reorder columns of x 
% giving x=[(x1(t-1) x2(t-1) .. xd(t-1)) (x1(t-2) x2(t-2)..xd(t-2)) ...
%           (x1(t-p) x2(t-p) .. xd(t-p))] on each row
for i=1:p
  start=(i-1)*d+1;
  stop=start+d-1;
  chunk(i).x=x(:,[start:1:stop]);
end
x=[];
for i=p:-1:1,
  x=[x chunk(i).x];
end
% and remove last row
Nrows=size(x,1);
x=x(1:Nrows-1,:);

y=X([p+1:1:N],:);

k=p*d*d;
if nargin < 3 || isempty(prior)
  prior.type='global';
  prior.groups=1;
  prior.group(1).index=ones(1,k);
end

% Get both pseudo-inverse and approx to inv(x'*x) efficiently
[ux,dx,vx]=svd(x);
kk=size(x,2);
ddx=diag(dx);
svd_tol=max(ddx)*eps*kk;
rank_X=sum(ddx > svd_tol);
ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';  % Pseudo-inverse
inv_xtx=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)'; % approx to inv(x'*x)

% Compute terms that will be used many times
xtx=x'*x;
yty=y'*y;
xty=x'*y;
vec_xty=xty(:);

% Get maximum likelihood solution
%w_ml = pinv(-1*x)*y;
w_ml = -xp*y;

% Swap signs to be consistent with paper (swap back at end !)
w_ml=-1*w_ml;

y_pred = x*w_ml;
e=y-y_pred;
noise_cov=(e'*e)/N;
sigma_ml=kron(noise_cov,inv_xtx);

% Priors on alpha(s)
b_alpha_prior=1000;
c_alpha_prior=0.001;
mean_alpha_prior=b_alpha_prior*c_alpha_prior;

% Initialise 
w_mean=w_ml;
w_cov=sigma_ml;

max_iters=32;
w=zeros(p*d,d);
tol=0.0001;
for it=1:max_iters,
    
    % Update weight precisions
    for j=1:prior.groups,
        Ij=diag(prior.group(j).index);
        kj=sum(prior.group(j).index);
        E_wj=0.5*w_mean(:)'*Ij*w_mean(:);
        b_alpha=E_wj+0.5*trace(Ij*w_cov*Ij)+(1/b_alpha_prior);
        group(j).b_alpha=1/b_alpha;
        group(j).c_alpha=0.5*kj+c_alpha_prior;
        group(j).mean_alpha=group(j).b_alpha*group(j).c_alpha;
        group(j).E_w=E_wj;
    end
    
    yy_pred=x*w_mean;
    ee=y-yy_pred;
    E_d_av=ee'*ee;
    
    Omega = spm_get_omega (p,d,w_cov,xtx);
    E_d_av=E_d_av+Omega;
    
    % Update noise precision posterior
    B=E_d_av;
    a=N;
    mean_lambda=a*inv(B);
    ilambda=(1/a)*B;
    
    prior_cov=zeros(k,k);
    for j=1:prior.groups,
        Ij=diag(prior.group(j).index);
        prior_cov=prior_cov+(1/group(j).mean_alpha)*Ij;
    end
    
    % Convergence criterion
    old_w=w;
    w=w_mean;
    
    if (it<=2)
        w=w_mean;
    else
        change=norm(w(:)-old_w(:))/k;
        if verbose
            fprintf('Iteration %d Delta_w=%1.6f',it,change);
        end
        if change < tol
            break;
        end
    end;
    
    % Update weight posterior
    data_precision=kron(mean_lambda,xtx);
    prior_prec=zeros(k,k);
    for j=1:prior.groups,
        Ij=diag(prior.group(j).index);
        prior_prec=prior_prec+group(j).mean_alpha*Ij;
    end
    w_cov=inv(data_precision+prior_prec);
    vec_w_mean=w_cov*data_precision*w_ml(:);
    w_mean=reshape(vec_w_mean,p*d,d);
    
end

% Compute Negative Free Energy

kl_alpha=0;
for j=1:prior.groups,
    kl_alpha=kl_alpha+spm_kl_gamma(group(j).b_alpha,group(j).c_alpha,b_alpha_prior,c_alpha_prior);
end
kl_weights=spm_kl_eig_normal(w_mean(:),w_cov,prior_cov);
lga=spm_lg_gamma(d,0.5*a);
acc=-0.5*N*log(det(B));
f_m=acc-kl_weights-kl_alpha+lga;
if verbose
    fprintf('Iteration %d Acc=%1.3f KL_w=%1.3f KL_alpha=%1.3f LogGena=%1.3f Fm=%1.3f',it,acc,kl_weights,kl_alpha,lga,f_m);
end

% Load up returning data structure
mar.p=p;

% This is the ML estimate
mar.noise_cov=noise_cov;

mar.sigma_ml=sigma_ml;

mar.w_cov=w_cov;
mar.group=group;
mar.prior=prior;
mar.mean_lambda=mean_lambda;
mar.fm=f_m;
mar.bic=-0.5*N*log(det(B))-0.5*k*log(N);
mar.iterations=it;
for i=1:p,
  start=(i-1)*d+1;
  stop=(i-1)*d+1+(d-1);
  % Transpose and swap signs for compatibility with spectral estimation function
  mar.lag(i).a=-w_mean(start:stop,:)';
end
% Swap signs for compatibility with spectral estimation function
mar.wmean=-w_mean;

% Compute the effective degrees of freedom
gamma=k;
for j=1:prior.groups,
    Ij=diag(prior.group(j).index);
    gamma=gamma-group(j).mean_alpha*trace(Ij*w_cov*Ij);
end
mar.gamma=gamma;
mar.d=d;

