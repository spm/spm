function [mar] = spm_csd2mar(csd,Hz,p,dt)
% Converts cross spectral density to MAR representation
% FORMAT [mar] = spm_csd2mar(csd,Hz,p,dt)
%
% csd  (N,:,:)   - cross spectral density
% Hz   (n x 1)   - vector of frequencies (Hz)
% p    (1)       - MAR(p) process    [default: p  = 8]
% dt             - sampling interval [default: dt = 1/(2*Hz(end))]
%
% mar  {1}       - see spm_mar
%
% See also: 
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m and spm_Q
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Defaults: MAR order and Nyquist
%--------------------------------------------------------------------------
if nargin < 3, p = 8;              end
if nargin < 4, dt = 1/(2*Hz(end)); end
 
% FFT and evaluate MAR coeficients
%--------------------------------------------------------------------------
ccf  = spm_csd2ccf(csd,Hz,dt);
mar  = spm_ccf2mar(ccf,p);
