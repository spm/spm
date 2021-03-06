function [y] = spm_fs_fmri_csd(y,M)
% Spectral feature selection for a CSD DCM
% FORMAT [y] = spm_fs_fmri_csd(y,M)
% y      - CSD
% M      - model structure
%__________________________________________________________________________
%
% This supplements cross spectral with cross covariance functions.
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% return (scaled) cross-spectra and covariance functions
%--------------------------------------------------------------------------
c    = spm_csd2ccf(y,M.Hz);
idx  = round(length(c(:,:,1))/2);
y    = [y; c(idx +(-8:8),:,:)*16];
