function [block] = spm_vb_alpha(Y,block)
% Variational Bayes for GLM-AR models - Update alpha 
% FORMAT [block] = spm_vb_alpha(Y,block)
%
% Y      - [T x N] time series 
% block  - data structure (see spm_vb_glmar)
%__________________________________________________________________________

% Will Penny and Nelson Trujillo-Barreto
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if block.verbose
    disp('Updating alpha');
end

N = block.N;
k = block.k;
  
% Convert from Sigma_n to Sigma_k
for n=1:N
    for j=1:k
        w_cov_k(n,j) = block.w_cov{n}(j,j);
    end
end
    
for j=1:k
    subblock_k = j:k:N*k;
    H = sum(spdiags(block.Dw,0).*w_cov_k(:,j)) + block.w_mean(subblock_k)'*block.Dw*block.w_mean(subblock_k);
    % Equation 15 in paper VB4
    block.b_alpha(j)    = 1./(H./2 + 1./block.b_alpha_prior(j));
    block.mean_alpha(j) = block.c_alpha(j)*block.b_alpha(j);
end
