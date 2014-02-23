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
% $Id: spm_fs_csd.m 5892 2014-02-23 11:00:16Z karl $


% return (scaled) cross-spectra and covariance functions
%--------------------------------------------------------------------------
p     = 16;
for i = 1:length(y);
    ccf  = spm_csd2ccf(y(i),M.Hz,M.dt);
    mar  = spm_ccf2mar(ccf{i},p);
    mar  = spm_mar_spectra(mar,M.Hz,1/M.dt);
    csd  = mar.P;
    ccf  = spm_csd2ccf(csd,M.Hz,M.dt);
    y{i} = [csd*8; ccf(1:8:end,:,:)];
end
