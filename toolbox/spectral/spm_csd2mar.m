function [mar] = spm_csd2mar(csd,Hz,p,dt)
% Converts cross spectral density to MAR representation
% FORMAT [mar] = spm_csd2mar(csd,Hz,p,dt)
%
% csd  (N,:,:)   - cross spectral density
% Hz   (n x 1)   - vector of frequencies (Hz)
% p    (1)       - MAR(p) process
% dt             - sampling interval
%
% amr  {1}       - see spm_mar
%
% See also: 
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m and spm_Q
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd2mar.m 5895 2014-02-26 14:28:23Z karl $

% Nyquist
%--------------------------------------------------------------------------
if nargin < 4, dt  = 1/(2*Hz(end)); end
 
% FFT and evalaute MAR coeficients
%--------------------------------------------------------------------------
ccf  = spm_csd2ccf(csd,Hz,dt);
mar  = spm_ccf2mar(ccf,p);




 
