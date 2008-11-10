function [block] = spm_vb_beta (Y,block)
% Variational Bayes for GLM-AR modelling in a block - Update beta 
% FORMAT [block] = spm_vb_beta (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_beta.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating beta');
end

N=block.N;
p=block.p;
  
% Convert from V_n to V_p
for n=1:N,
    for j = 1:p,
        a_cov_p(n,j) = block.a_cov{n}(j,j);
    end
end
    
switch block.priors.A,
    case 'Discrete',
        for j = 1:p,
            for s=1:block.priors.S,
                sj=(block.priors.voxel(s).i-1)*p+j;
                H=sum((block.a_mean(sj)-block.as(j,s)).^2+a_cov_p(block.priors.voxel(s).i,j));
                block.b_beta(j,s)=1/(H/2+1./block.b_beta_prior(j,s));
                block.mean_beta(j,s) = block.c_beta(j,s)*block.b_beta(j,s);
            end
        end
    otherwise
        for j = 1:p,
            subblock_p = j:p:N*p;
            H  = sum(spdiags(block.Da,0).*a_cov_p(:,j)) + block.a_mean(subblock_p)'*block.Da*block.a_mean(subblock_p);
            % Equation 16 in paper VB4
            block.b_beta(j)    = 1./(H./2 + 1./block.b_beta_prior(j));
            block.mean_beta(j) = block.c_beta(j)*block.b_beta(j);
        end
end