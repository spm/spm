function [ccf,pst] = spm_ker2ccf(ker,dt)
% computes cross covariance function from kernels
% FORMAT [ccf,pst] = spm_ker2ccf(ker,dt)
%
% ker  - first-order (Volterra) kernels
% dt   - time bin (sec)
%
% ccf  - cross covariance functions
% pst  - time samples
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% cross covariance function
%==========================================================================

% convolution of kernels
%--------------------------------------------------------------------------
N     = 2*size(ker,1) - 1;
pst   = (1:N) - mean(1:N);
pst   = pst*dt;
ccf   = zeros(N,size(ker,2),size(ker,2));
for i = 1:size(ker,2)
    for j = 1:size(ker,2)
        for k = 1:size(ker,3)
            ccf(:,i,j) = ccf(:,i,j) + conv(ker(:,i,k),flip(ker(:,j,k)));
        end
    end
end
ccf   = ccf*dt/2;


% via modulation transfer function
%--------------------------------------------------------------------------
% [mtf,Hz]  = spm_ker2mtf(ker,dt);
% [ccf,pst] = spm_mtf2ccf(mtf,Hz);