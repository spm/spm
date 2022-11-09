function [csd,Hz] = spm_ker2csd(ker,pst)
% computes cross spectral density from kernels
% FORMAT [csd,Hz] = spm_ker2csd(ker,pst)
%
% ker  - first-order (Volterra) kernels
% pst  - time samples
%
% csd  - cross spectral density
% Hz   - frequencies
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% cross spectral density
%==========================================================================

% via modulation transfer function
%--------------------------------------------------------------------------
[mtf,Hz] = spm_ker2mtf(ker,pst);
csd      = spm_mtf2csd(mtf);