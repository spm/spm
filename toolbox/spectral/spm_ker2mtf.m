function [mtf,Hz] = spm_ker2mtf(ker,dt)
% computes modulation transfer function from kernels
% FORMAT [mtf,Hz] = spm_ker2mtf(ker,dt)
%
% ker  - first-order (Volterra) kernels
% dt   - time bin
%
% mtf  - modulation transfer function
% Hz   - frequencies
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------


% Transfer functions
%==========================================================================

% Transfer functions (FFT of kernel)
%--------------------------------------------------------------------------
N    = size(ker,1);
mtf  = fft(ker,2*N + 1)*dt;
Hz   = ((1:N) - 1)/(N*dt)/2;
