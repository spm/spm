function [csd,Hz] = spm_ker2csd(ker,pst)
% computes cross spectral density from kernels
% FORMAT [csd,Hz] = spm_ker2csd(ker,pst)
%
% ker  - first-order (Volterra) kernels
% pst  - time samples
%
% csd  - cross spectral density
% Hz   - frequencies
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ker2csd.m 7774 2020-01-25 18:07:03Z karl $


% cross spectral density
%==========================================================================

% via modulation transfer function
%--------------------------------------------------------------------------
[mtf,Hz] = spm_ker2mtf(ker,pst);
csd      = spm_mtf2csd(mtf);