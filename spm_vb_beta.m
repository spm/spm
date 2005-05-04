function [slice] = spm_vb_beta (Y,slice)
% Variational Bayes for GLM-AR modelling in a slice - Update beta 
% FORMAT [slice] = spm_vb_beta (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure 
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_beta.m 112 2005-05-04 18:20:52Z john $

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
    
switch slice.priors.A,
    case 'Discrete',
        for j = 1:p,
            for s=1:slice.priors.S,
                sj=(slice.priors.voxel(s).i-1)*p+j;
                H=sum((slice.a_mean(sj)-slice.as(j,s)).^2+a_cov_p(slice.priors.voxel(s).i,j));
                slice.b_beta(j,s)=1/(H/2+1./slice.b_beta_prior(j,s));
                slice.mean_beta(j,s) = slice.c_beta(j,s)*slice.b_beta(j,s);
            end
        end
    otherwise
        for j = 1:p,
            block_p = j:p:N*p;
            H  = sum(spdiags(slice.Da,0).*a_cov_p(:,j)) + slice.a_mean(block_p)'*slice.Da*slice.a_mean(block_p);
            % Equation 16 in paper VB4
            slice.b_beta(j)    = 1./(H./2 + 1./slice.b_beta_prior(j));
            slice.mean_beta(j) = slice.c_beta(j)*slice.b_beta(j);
        end
end