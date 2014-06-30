function [block] = spm_vb_w(Y,block)
% Variational Bayes for GLM-AR modelling in a block - update w
% FORMAT [block] = spm_vb_w (Y,block)
%
% Y          - [T x N] time series 
% block      - data structure (see spm_vb_glmar)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_w.m 6079 2014-06-30 18:25:37Z spm $

if block.verbose
    disp('Updating w');
end

X = block.X;
T = block.T;
p = block.p;
N = block.N;
k = block.k;
Bk = kron(diag(block.mean_alpha),block.Dw);
B = block.Hw*Bk*block.Hw';

if p > 0
    voxel = spm_vb_get_Ab (Y,block);
end

for n = 1:N,
    subblock_n              = [(n-1)*k+1:n*k];
    subblock_ni             = [1:N*k];
    subblock_ni(subblock_n) = [];
    Bnn                     = B(subblock_n,subblock_n);
    Bni                     = B(subblock_n,subblock_ni);
        
    if p > 0
        block.w_cov{n}  = inv(block.mean_lambda(n)*voxel(n).A + Bnn);
        w_mean = block.w_cov{n} * (block.mean_lambda(n)*voxel(n).b-Bni*block.w_mean(subblock_ni,1));
    else
        block.w_cov{n}  = inv(block.mean_lambda(n)*block.XTX + Bnn);
        w_mean = block.w_cov{n} * (block.mean_lambda(n)*block.XTY(:,n)-Bni*block.w_mean(subblock_ni,1));
    end
    block.w_mean(subblock_n,1) = w_mean;
    block.w2{n} = w_mean*w_mean' + block.w_cov{n};
    
    % Update intermediate quantities
    if p > 0
        block.I.W_tilde(:,:,n) = reshape(block.I.D(:,:,n)'*w_mean,p,p);
    end
end
