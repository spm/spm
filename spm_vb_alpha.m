function [slice] = spm_vb_alpha (Y,slice)
% Variational Bayes for GLM-AR models - Update alpha 
% FORMAT [slice] = spm_vb_alpha (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure 
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_alpha.m 112 2005-05-04 18:20:52Z john $

if slice.verbose
    disp('Updating alpha');
end

N=slice.N;
k=slice.k;
  
% Convert from Sigma_n to Sigma_k
for n=1:N,
    for j = 1:k,
        w_cov_k(n,j) = slice.w_cov{n}(j,j);
    end
end
    
for j = 1:k,
    block_k = j:k:N*k;
    H  = sum(spdiags(slice.Dw,0).*w_cov_k(:,j)) + slice.w_mean(block_k)'*slice.Dw*slice.w_mean(block_k);
    % Equation 15 in paper VB4
    slice.b_alpha(j)    = 1./(H./2 + 1./slice.b_alpha_prior(j));
    slice.mean_alpha(j) = slice.c_alpha(j)*slice.b_alpha(j);
end
         
    