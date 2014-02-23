function [gew] = spm_ccf2gew(ccf,Hz,dt,p)
% Converts cross covariance function to Geweke Granger causality
% FORMAT [gew] = spm_ccf2gew(ccf,Hz,dt,p)
%
% ccf  (N,m,m)   - cross covariance functions
% Hz   (n x 1)   - vector of frequencies (Hz)
% dt             - samping interval
% p              - AR(p) order [default p = 8]
%
% gwe   (N,m,m)  - Geweke's frequency domain Granger causality
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ccf2gew.m 5892 2014-02-23 11:00:16Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 4, p = 8; end

% Granger causality
%==========================================================================
mar  = spm_ccf2mar(ccf,p);
mar  = spm_mar_spectra(mar,Hz,1/dt);
gew  = mar.gew;