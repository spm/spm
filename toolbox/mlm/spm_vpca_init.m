function [W_ml,lambda,sigma2] = spm_vpca_init (T, form_cov)
% Initialise VPCA model
% function [W_ml,lambda,sigma2] = spm_vpca_init (T, form_cov)
%
% T         [d x N] matrix containing N d-dimensional data vectors
%           The nth data sample, t_n, is nth column of T
%
% form_cov  form covariance matrix (1=yes, 0=no, default=no)
% 
% W_ml      Maximum Likelihood (ML) estimate of factor matrix
% lambda    eigenvalues
% sigma2    Observation noise variance
%__________________________________________________________________________
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vpca_init.m 5962 2014-04-17 12:47:43Z spm $

% Maximum size of data set that SVD can handle given typical PC memory
max_size=10^6;

[d,N]=size(T);

if d*N > max_size
    % Random initialisation
    q=min(d,N)-1;
    W_ml=randn(d,q);
    lambda=1;
    sigma2=1;
    return
end

q=min(d,N)-1; % Set latent space dimensionality to max possible 
noise_dim=q+1;

if nargin < 2 || isempty(form_cov)
    form_cov=0;
end

if form_cov
    S=cov(T');
    [v,ds]=eig(S);
    dd=diag(ds);
    [ddd,i]=sort(dd);
    i=flipud(i);
    v=v(:,i);
    lambda=dd(i);
    v=v(:,1:q);    
else
    % Economy size SVD decomposition
    [U,S,V]=svd(sqrt(1/N)*T,0);
    v=U(:,1:q);
    s=diag(S);
    lambda=s.^2;
end
sigma2=mean(lambda(noise_dim));
lambda=lambda(1:q);
    
W_ml=v*diag(sqrt(lambda-sigma2));
