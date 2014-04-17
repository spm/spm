function [pca] = spm_vpca (T,q,Bayes)
% Variational Principal Component Analysis 
% FORMAT [pca] = spm_vpca (T,q,Bayes)
%
% T     [d x N] matrix containing N d-dimensional data vectors
%       The nth data sample, t_n, is nth column of T
% q     maximum latent space dimension (q < d)
% Bayes 1 for Bayesian algorithm, 0 otherwise (default = 1)
%
% pca model is
%
%       t_n = W x_n + mu + e
%
% See C. Bishop. Variational Principal Components, ANN, 1999.
%
% The factor matrix W is a [d x q] matrix, where q=d-1
% The ith factor is in ith column
%
% pca  Contains fields for
%
%       ML solution: ml.W, ml.lambda (factor matrix and eigenvalues)
%       Latent variables: M_x, Sigma_x
%       Mean: mean_mu, Sigma_mu
%       Factor Matrix: M_w, Sigma_w 
%       Predicted Data: That, mse (mean square error of predictions)
%       Neg Free Energy: Fm, Fm_evol
%       Observation noise precision: mean_tau
%       Prior precisions of Factor magnitudes: mean_alpha
%       Prior precision of Mean: mean_beta
%
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vpca.m 5962 2014-04-17 12:47:43Z spm $

[d,N]=size(T);
if nargin < 2 || isempty(q)
    % Set latent space dimensionality to max possible 
    q=d-1;
end
if nargin < 3 || isempty(Bayes)
    Bayes=1;
end

% Set priors
a_alpha=0.001;
b_alpha=0.001;
a_tau=0.001;
b_tau=0.001;
mean_tau=a_tau/b_tau;

% Initialise component
c(1).mean_mu=mean(T,2);
[W_ml,lambda] = spm_vpca_init (T-c(1).mean_mu*ones(1,N));
c(1).q=q;
c(1).M_w=W_ml(:,1:q);
c(1).avg_WtW=c(1).M_w'*c(1).M_w;
c(1).mean_alpha=(std(c(1).M_w,[],1)+eps).^(-2);
        
pca.ml.W=W_ml;
pca.ml.lambda=lambda;
if ~Bayes
    return
end

% prior precision of means - may need to adjust to scale of data
beta=0.001; 
ibeta=1/beta;

S=ones(1,N);
    
% Initialise
tic;
pca.N=N;pca.d=d;
pca.a_alpha=a_alpha;
pca.b_alpha=b_alpha;
pca.beta=beta;
pca.a_tau=a_tau;
pca.b_tau=b_tau;
pca.mean_tau=mean_tau;

pca.lambda=N;
pca.u=1;
pca.M=1;
Nm=N;
    
pca_T=T;
pca_S=S;
pca_c=c;

maxits=32;
Fm_evol=[];
for it=1:maxits,
    
    [pca,c] = spm_vpca_update (T,S,pca,c,1);
    
    Fm = spm_vpca_f (pca,c);
        
    if it > 1
        Fm_evol=[Fm_evol,Fm];
    end
    if it > 1
        disp(sprintf('PCA iteration %d: F = %1.2f',it,Fm));
    end
    
end % Iteration of update equations

% Get MSE of model predictions
for n=1:N,
    mse(:,n)=zeros(d,1);
    That(:,n)=zeros(d,1);
    Tpred=c(1).M_w*c(1).M_x(:,n)+c(1).mean_mu;
    mnerr=(T(:,n)-Tpred);
    mse(:,n)=mse(:,n)+mnerr.^2;
    That(:,n)=That(:,n)+Tpred;
end
mse=mean(mean(mse));

% Fill in o/p data structure
pca.Fm=Fm;
pca.Fm_evol=Fm_evol;

pca.M_x=c(1).M_x;
pca.Sigma_x=c(1).Sigma_x;
pca.mean_mu=c(1).mean_mu;
pca.Sigma_mu=c(1).Sigma_mu;
pca.M_w=c(1).M_w;
pca.Sigma_w=c(1).Sigma_w;
pca.mean_alpha=c(1).mean_alpha;
pca.mean_tau=pca.qa_tau/pca.qb_tau;
pca.mse=mse;
pca.That=That;
