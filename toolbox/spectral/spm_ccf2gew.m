function [gew] = spm_ccf2gew(ccf,Hz,dt,p)
% Converts cross covariance function to Geweke Granger causality
% FORMAT [gew] = spm_ccf2gew(ccf,Hz,dt,p)
%
% ccf  (N,m,m)   - cross covariance functions
% Hz   (n x 1)   - vector of frequencies (Hz)
% dt             - samping interval [default dt = 1/(2*Hz(end))]
% p              - AR(p) order [default p = 8]
%
% gwe   (N,m,m)  - Geweke's frequency domain Granger causality
%
% See also: spm_???2???.m
%     ??? = {'ccf','csd','gew','mar','coh','mtf','ker','ssm','dcm'}
% and spm_Q.m, spm_mar.m, spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ccf2gew.m 7774 2020-01-25 18:07:03Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3, dt = 1/(2*Hz(end)); end                % Nyquist
if nargin < 4, p  = 8;             end                % MAR order

% Granger causality
%==========================================================================
mar  = spm_ccf2mar(ccf,p);
mar  = spm_mar_spectra(mar,Hz,1/dt);
gew  = mar.gew;