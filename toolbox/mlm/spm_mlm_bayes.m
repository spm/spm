function [mlm] = spm_mlm_bayes (y,x,pr,verbose,ml_only)
% Bayesian Multivariate Linear Modelling
% FORMAT [mlm] = spm_mlm_bayes (y,x,pr,verbose,ml_only)
%
% MLM: y = x W + e
%
% y           T-by-d data matrix 
% x           N-by-p design matrix
% pr          Shrinkage prior on MLM coefficients:
%             'input' (default), 'output' or 'global'
%
%             For 'input', coeffs of each independent variable
%             ie rows of W, share same prior precision. This 
%             allows some inputs to be more relevant than others.
%
%             For 'output', cols of W share same prior precision.
%             This allows some outputs to be more relevant.
%
%             For 'global' there is a single prior precision.
%
% verbose     1 to print out iteration details, 0 otherwise (default=0)
% ml_only     set to 1 to only compute ML solution. Default is zero
%
% The returned data structure mlm contains the following fields
%
% .wmean      Bayes estimate of [p x d] regression coefficient matrix
% .wsd        [p x d] posterior standard deviations of reg coeffs
% .wml        Maximum Likelihood regression coefficient matrix
% .wcov       [pd x pd] posterior covariance of regression coeffs
% .lambda     [d x d] observation noise precision matrix
% .fm         Negative free energy of model
% .bic        Bayesian Information Criterion
% .iterations Number of iterations during optimisation
% .prior      Details of regression coeff prior
%             .group(j).mean_alpha:
%             Estimated prior precision of jth parameter group.
%             For 'input' prior this is jth row of W. 
%             For 'output' prior this is jth column of W.
%___________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mlm_bayes.m 4651 2012-02-09 16:03:39Z will $

if nargin < 4 || isempty(verbose)
    verbose=0;
end
if nargin < 5 || isempty(ml_only)
    ml_only=0;
end

d=size(y,2);    % Dimension of each dependent variable
N=size(y,1);    % Number of data points
p=size(x,2);    % Dimension of independent variable

if d >= N | p >= N
    disp('Error in spm_mlm_bayes.m: too few data points');
    return
end

k=p*d;
if nargin < 3 || isempty(pr)
    pr='input';
end

switch lower(pr)
    case 'input',
        % Separate precision for each row of W ie for each input
        prior.type='input';
        prior.groups=p;
        % get indices of each 'parameter group' ie rows
        vec_ind=[1:k]';
        ind=reshape(vec_ind,p,d);
        g=zeros(1,k);
        for j=1:p,
            gj=g; gj(ind(j,:))=1;
            group(j).index=gj;
        end
        
    case 'output',
        % Separate precision for each column of W ie for each output
        prior.type='output';
        prior.groups=d;
        % get indices of each 'parameter group' ie cols
        vec_ind=[1:k]';
        ind=reshape(vec_ind,p,d);
        g=zeros(1,k);
        for j=1:d,
            gj=g; gj(ind(:,j))=1;
            group(j).index=gj;
        end
        
    otherwise
        % Global prior
        prior.type='global';
        prior.groups=1;
        group(1).index=ones(1,k);
end

% Get both pseudo-inverse and approx to inv(x'*x) efficiently
[ux,dx,vx]=svd(x);
kk=size(x,2);
if kk==1
    xp=pinv(x);
    inv_xtx=1/(x'*x);
else
    ddx=diag(dx);
    svd_tol=max(ddx)*eps*kk;
    rank_X=sum(ddx > svd_tol);
    ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
    ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
    xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';  % Pseudo-inverse
    inv_xtx=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)'; % approx to inv(x'*x)
end

% Compute terms that will be used many times
xtx=x'*x;
yty=y'*y;
xty=x'*y;
vec_xty=xty(:);

% Get maximum likelihood solution
w_ml = xp*y;
if ml_only
    mlm.wml=w_ml;
    return
end

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
w=zeros(p,d);
tol=0.0001;
for it=1:max_iters,
    
    % Update weight precisions
    for j=1:prior.groups,
        Ij=diag(group(j).index);
        kj=sum(group(j).index);
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
    
    Omega = get_omega (p,d,w_cov,xtx);
    E_d_av=E_d_av+Omega;
       
    % Update noise precision posterior
    B=E_d_av;
    a=N;
    mean_lambda=a*inv(B);
    ilambda=(1/a)*B;
    
    prior_cov=zeros(k,k);
    for j=1:prior.groups,
        Ij=diag(group(j).index);
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
            disp(sprintf('Iteration %d Delta_w=%1.6f',it,change));
        end
        if change < tol
            break;
        end
    end;
   
    % Update weight posterior
    data_precision=kron(mean_lambda,xtx);
    prior_prec=zeros(k,k);
    for j=1:prior.groups,
        Ij=diag(group(j).index);
        prior_prec=prior_prec+group(j).mean_alpha*Ij;
    end
    w_cov=inv(data_precision+prior_prec);
    vec_w_mean=w_cov*data_precision*w_ml(:);
    
    w_mean=reshape(vec_w_mean,p,d);
    
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
    disp('Contributions to Negative Free Energy');
    disp(sprintf('Acc=%1.3f KL_w=%1.3f KL_alpha=%1.3f LogGena=%1.3f Fm=%1.3f',acc,kl_weights,kl_alpha,lga,f_m));
end

% Get error bars on regression coefficients
w_sd=zeros(p,d);
post_var=diag(w_cov);
for dd=1:d,
    start=(dd-1)*p+1;
    stop=start+p-1;
    w_sd(:,dd)=sqrt(post_var(start:stop));
end

% Create mlm data structure
mlm.wmean=w_mean;
mlm.wsd=w_sd;
mlm.wml=w_ml;
mlm.wcov=w_cov;
mlm.lambda=mean_lambda;
mlm.fm=f_m;
mlm.bic=-0.5*N*log(det(B))-0.5*k*log(N);
mlm.iterations=it;

mlm.prior=prior;
mlm.prior.group=group;


function [Omega] = get_omega (p,d,w_cov,xtx)
% Get contribution to prediction error variance from w_cov 
% FORMAT [Omega] = spm_mlm_get_omega (p,d,w_cov,xtx)
%
% p         Number of independent variables
% d         Number of dependent variables
% w_cov     Uncertainty in MLM coefficients
% xtx       X'X where X is design matrix 
%
% Omega     Expected error variance from w_cov


Omega=zeros(d,d);
% Submatrix size - ie number of model params per dependent variable
s=p;

% Get upper diagonal elements
for di=1:d,
    for dj=di:d,
        istart=1+(di-1)*s;
        istop=istart+s-1;
        jstart=1+(dj-1)*s;
        jstop=jstart+s-1;
        w_cov_i_j=w_cov(istart:istop,jstart:jstop);
        Omega(di,dj)=trace(w_cov_i_j*xtx);
    end
end

% Get lower diagonal elements
Omega = Omega+Omega'-diag(diag(Omega));