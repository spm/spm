function [ccf,pst] = spm_ker2ccf(ker,dt)
% computes cross covariance function from kernels
% FORMAT [ccf,pst] = spm_ker2ccf(ker,dt)
%
% ker  - first-order (Volterra) kernels
% dt   - time bin (sec)
%
% ccf  - cross covariance functions
% pst  - time samples
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ker2ccf.m 7774 2020-01-25 18:07:03Z karl $


% cross covariance function
%==========================================================================

% via modulation transfer function
%--------------------------------------------------------------------------
[mtf,Hz]  = spm_ker2mtf(ker,dt);
[ccf,pst] = spm_mtf2ccf(mtf,Hz);