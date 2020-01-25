function [mtf,Hz] = spm_ker2mtf(ker,dt)
% computes modulation transfer function from kernels
% FORMAT [mtf,Hz] = spm_ker2mtf(ker,dt)
%
% ker  - first-order (Volterra) kernels
% dt   - time bin
%
% mtf  - modulation transfer function
% Hz   - frequencies
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ker2mtf.m 7774 2020-01-25 18:07:03Z karl $

% preliminaries
%--------------------------------------------------------------------------


% Transfer functions
%==========================================================================

% Transfer functions (FFT of kernel)
%--------------------------------------------------------------------------
N    = size(ker,1);
mtf  = fft(ker,2*N + 1)*dt;
Hz   = ((1:N) - 1)/(N*dt)/2;





