function [y] = spm_fs_csd(y,M)
% spectral feature selection for a CSD DCM
% FORMAT [y] = spm_fs_csd(y,M)
% y      - CSD
% M      - model structure
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fs_csd.m 6557 2015-09-20 12:44:30Z karl $


% return (scaled) cross-spectra and covariance functions
%--------------------------------------------------------------------------
for i = 1:length(y);
    csd  = y{i};
    ccf  = spm_csd2ccf(csd,M.Hz);
    y{i} = [csd*8; ccf(1:8:end,:,:)];
    % y{i} = [log(csd); ccf(1:8:end,:,:)];
end
