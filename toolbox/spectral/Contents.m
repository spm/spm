% SPM Spectral Estimation Toolbox
%
% Installation note: Some of these routines require the SPM Mixture 
% Modelling Toolbox to be on the search path
%
% Bayesian Autoregressive (AR) modelling: spm_ar.m
% 
% Bayesian Robust Autoregressive (AR) modelling: spm_rar.m
%
% Bayesian Multivariate Autoregressive (MAR) modelling [1,2,3]: spm_mar.m
%
% (MAR based) Granger causality [3]: spm_granger.m
%
% (MAR based) Spectral estimation: spm_mar_spectra.m
%
% Wavelet based spectrogram: spm_wavspec.m
%
% References:
% 
% [1] W.D. Penny and S.J. Roberts. Bayesian Multivariate Autoregresive Models 
% with structured priors. IEE Proceedings on Vision, Image and Signal Processing, 149(1):33-41, 2002
%
% [2] L. Harrison, W.D. Penny, and K.J. Friston. Multivariate Autoregressive 
% Modelling of fMRI time series. NeuroImage, 19(4):1477-1491, 2003
% 
% [3] W. Penny and L. Harrison. Multivariate autoregressive models. In K. Friston, 
% J. Ashburner, S. Kiebel, T. Nichols, and W. Penny, editors, Statistical 
% Parametric Mapping: The analysis of functional brain images. Elsevier, 
% London, 2006

%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: Contents.m 1143 2008-02-07 19:33:33Z spm $