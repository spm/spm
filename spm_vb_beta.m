function [slice] = spm_vb_beta (Y,slice)
% Variational Bayes for GLM-AR modelling in a slice - Update beta 
% FORMAT [slice] = spm_vb_beta (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure 
%
% %W% Will Penny and Nelson Trujillo-Barreto %E%

if slice.verbose
    disp('Updating beta');
end

N=slice.N;
p=slice.p;
  
% Convert from V_n to V_p
for n=1:N,
    for j = 1:p,
        a_cov_p(n,j) = slice.a_cov{n}(j,j);
    end
end
    
for j = 1:p,
    block_p = j:p:N*p;
    H  = sum(spdiags(slice.D,0).*a_cov_p(:,j)) + slice.a_mean(block_p)'*slice.D*slice.a_mean(block_p);
    % Equation 16 in paper VB4
    slice.b_beta(j)    = 1./(H./2 + 1./slice.b_beta_prior(j));
    slice.mean_beta(j) = slice.c_beta(j)*slice.b_beta(j);
end
         
    