function [M_opt,log_ev,lambda,var] = spm_pca_order (X, N)
% Model order selection for PCA   
% FORMAT [M_opt,log_ev,lambda,var] = spm_pca_order (X, N)
%
% Model order selection for PCA using Minka's approximation to model evidence
% Input can be
%     X         Data
% or
%    
%     X         Covariance matrix
%     N         number of samples used for computing X
%
% M_opt         Optimum number of sources
% log_ev        Log Evidence
% lambda        Eigenspectrum
% var           Estimated observation noise (at M_opt)
%
% Algorithm:
%
% T.P. Minka. Automatic choice of dimensionality for PCA. Technical Report
% 514, MIT Media Lab, Perceptual Computing Section, 2000.
%
% Evaluation:
%
% W. Penny, S. Roberts and R. Everson (2000) ICA: model order selection
% and dynamic source models. ICA: Principles and Practice, pages 299-314. 
% Cambridge University Press.
%__________________________________________________________________________
% Copyright (C) 2007-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_pca_order.m 5962 2014-04-17 12:47:43Z spm $


if nargin == 1
    [N,d]=size(X);
    if d > N
        X=X';
        [N,d]=size(X);
    end
    X=X-ones(N,1)*mean(X);
    S=(1/N)*ctranspose(X)*X; %to also work for complex input
else
    d = size(X, 1);
    S = X;
end

[w,lambda] = eig (S);
lambda=diag(lambda);

% Order eigenvectors/values
lambda=sort(lambda);
lambda=flipud(lambda);

% Loop over possible number of sources
for M=1:d-1,
  % Minka equation 50
  i=[1:1:M];
  kk=(d-i+1)/2;
  term1=-M*log(2)+sum(gammaln(kk))+sum(-kk*log(pi));
  term2=-0.5*N*sum(log(lambda(i)));
  var=mean(lambda(M+1:d));
  term3=-0.5*N*(d-M)*log(var);
  little_m=d*M-M*(M+1)/2;
  term4=0.5*(little_m+M)*log(2*pi);
  
  lambda_hat=[lambda(1:M);var*ones(d-M,1)];
  term5=0;
  for i=1:M,
    for j=i+1:d,
      term5=term5+log(1/lambda_hat(j)-1/lambda_hat(i))+log(lambda(i)-lambda(j))+ log(N);
    end
  end
  term5=-0.5*term5;
  term6=-0.5*M*log(N);

  % Minka equation 73
  log_ev(M)=term1+term2+term3+term4+term5+term6;
end

[max_ev,M_opt]=max(log_ev);

var=mean(lambda(M_opt+1:d));
