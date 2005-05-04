function [slice] = spm_vb_gamma (Y,slice)
% Variational Bayes for GLMAR model - Update gamma
% and get w_dev, wk_mean
% FORMAT [slice] = spm_vb_gamma (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure containing the following fields:
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_gamma.m 112 2005-05-04 18:20:52Z john $

if slice.verbose
    disp('Updating gamma');
end

N = slice.N;
k = slice.k;
k = slice.k;

Bk = kron(diag(slice.mean_alpha),slice.Dw);
B = slice.Hw*Bk*slice.Hw';

for n=1:N,
    % Block matrices Bnn [k x k] and Bni [k x k*(N-1)]
    block_n           = [(n-1)*k+1:n*k];
    Bnn               = B(block_n,block_n);
    % Equation 17 in paper VB2
    for j=1:k,
        slice.gamma(j,n)= 1-slice.w_cov{n}(j,j)*Bnn(j,j);
        slice.b(j,n)=Bnn(j,j);
    end
    % Record Standard Deviation of parameter estimates
    % to be used in Taylor series approximation to posterior
    slice.w_dev(:,n) = sqrt(diag(slice.w_cov{n}));
end
slice.gamma_tot=sum(slice.gamma,2);
slice.wk_mean      = reshape(slice.w_mean,k,N); 
