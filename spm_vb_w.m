function [slice] = spm_vb_w (Y,slice)
% Variational Bayes for GLM-AR modelling in a slice - update w
% FORMAT [slice] = spm_vb_w (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure 
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_w.m 112 2005-05-04 18:20:52Z john $

if slice.verbose
    disp('Updating w');
end

X = slice.X;
T = slice.T;
p = slice.p;
N = slice.N;
k = slice.k;
Bk = kron(diag(slice.mean_alpha),slice.Dw);
B = slice.Hw*Bk*slice.Hw';

if p > 0
    voxel=spm_vb_get_Ab (Y,slice);
end

for n=1:N,
    block_n           = [(n-1)*k+1:n*k];
    block_ni          = [1:N*k];
    block_ni(block_n) = [];
    Bnn               = B(block_n,block_n);
    Bni               = B(block_n,block_ni);
        
    if p > 0
        slice.w_cov{n}        = inv(slice.mean_lambda(n)*voxel(n).A + Bnn);
        w_mean = slice.w_cov{n} * (slice.mean_lambda(n)*voxel(n).b-Bni*slice.w_mean(block_ni,1));
    else
        slice.w_cov{n}        = inv(slice.mean_lambda(n)*slice.XTX + Bnn);
        w_mean = slice.w_cov{n} * (slice.mean_lambda(n)*slice.XTY(:,n)-Bni*slice.w_mean(block_ni,1));
    end
    slice.w_mean(block_n,1)=w_mean;
    slice.w2{n}=w_mean*w_mean'+slice.w_cov{n};
    
    % Update intermediate quantities
    if p>0
        slice.I.W_tilde(:,:,n)=reshape(slice.I.D(:,:,n)'*w_mean,p,p);
    end
end
    
