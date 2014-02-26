function [y] = spm_fs_csd(y,M)
% spectral feature selection for a CSD DCM
% FORMAT [y] = spm_fs_csd(y,M)
% y      - CSD
% M      - model structure
%__________________________________________________________________________
%
% This simply log-transforms the (real) auto-spectra
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fs_csd.m 5895 2014-02-26 14:28:23Z karl $


% return (scaled) cross-spectra and covariance functions
%--------------------------------------------------------------------------
for i = 1:length(y);
    ccf  = spm_csd2ccf(y{i},M.Hz);
    y{i} = [y{i}*8; ccf(1:8:end,:,:)];
end
