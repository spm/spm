function [block] = spm_vb_gamma(Y,block)
% Variational Bayes for GLMAR model - Update gamma and get w_dev, wk_mean
% FORMAT [block] = spm_vb_gamma(Y,block)
%
% Y      - [T x N] time series
% block  - data structure (see spm_vb_glmar)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_gamma.m 6079 2014-06-30 18:25:37Z spm $

if block.verbose
    disp('Updating gamma');
end

N  = block.N;
k  = block.k;
k  = block.k;

Bk = kron(diag(block.mean_alpha),block.Dw);
B  = block.Hw*Bk*block.Hw';

for n=1:N
    % Block matrices Bnn [k x k] and Bni [k x k*(N-1)]
    subblock_n           = [(n-1)*k+1:n*k];
    Bnn                  = B(subblock_n,subblock_n);
    % Equation 17 in paper VB2
    for j=1:k
        block.gamma(j,n) = 1-block.w_cov{n}(j,j)*Bnn(j,j);
        block.b(j,n)     = Bnn(j,j);
    end
    % Record Standard Deviation of parameter estimates
    % to be used in Taylor series approximation to posterior
    block.w_dev(:,n)     = sqrt(diag(block.w_cov{n}));
end
block.gamma_tot          = sum(block.gamma,2);
block.wk_mean            = reshape(block.w_mean,k,N); 
