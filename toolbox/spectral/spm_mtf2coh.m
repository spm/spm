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

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% convert via cross spectral density
%==========================================================================
csd       = spm_mtf2csd(mtf);
[coh,fsd] = spm_csd2coh(csd,Hz);
