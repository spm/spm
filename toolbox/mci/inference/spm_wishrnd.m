function [S] = spm_wishrnd (B,a,N)
% Generate N samples from Wishart density
% FORMAT [S] = spm_wishrnd (B,a,N)
%
% B,a   Wishart params, d=dim(B)
% N     Number of samples
% S     [d x d x N] sample matrices or [d x d] if N=1
%
% The Wishart density here, W(S;a,B), is defined as in p. 435 of
% J. Bernardo and A. Smith, Bayesian Theory, Wiley, 2000. 
% We have E[S]=aB^{-1}
% 
% This definition is different to eg. C. Bishop, 
% Pattern Recognition and Machine Learning, Springer, 2006., who
% have W(S;n,V). They are related by n=2a, V=B^{-1}/2. We have E[S]=nV
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_wishrnd.m 6548 2015-09-11 12:39:47Z will $

if nargin < 3 | isempty(N)
    N=1;
end

d = size(B,1);
n = 2*a;
V = 0.5*inv(B);

% Generate scatter matrices i.e. outer product of n Gaussian variates 
% with mean 0 and covariance V
x = spm_normrnd(zeros(d,1),V,n*N);
S = zeros(d,d,N);
for s=1:N,
    ind = [(s-1)*n+1:s*n];
    S(:,:,s) = x(:,ind)*x(:,ind)';
end

if N==1
    S = squeeze(S);
end
        