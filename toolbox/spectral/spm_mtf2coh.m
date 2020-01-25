function [coh,fsd] = spm_mtf2coh(mtf,Hz)
% Converts modulation transfer function to coherence
% FORMAT [coh,fsd] = spm_mtf2coh(mtf,Hz)
%
% mtf  (N,n,n)   - (unnormalised) directed or modulation transfer function
% Hz   (N x 1)   - vector of frequencies (Hz)
%
% coh            - coherence
% fsd            - frequency specific delay (seconds) 
%                - phase-delay/radial frequency
%
% See also: 
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mtf2coh.m 7774 2020-01-25 18:07:03Z karl $

% convert via cross spectral density
%==========================================================================
csd       = spm_mtf2csd(mtf);
[coh,fsd] = spm_csd2coh(csd,Hz);