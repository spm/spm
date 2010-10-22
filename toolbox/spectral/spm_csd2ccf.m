function [ccf,pst] = spm_csd2ccf(csd,Hz)
% Converts cross spectral density to cross covariance function
% FORMAT [ccf,pst] = spm_csd2ccf(csd,Hz)
%
% csd  (Hz,:,:)         - cross spectral density (cf, mar.P)
% Hz   (n x 1)          - vector of frequencies (Hertz)
%
% ccf                   - cross covariance functions
% pst  (n + n + 1 x 1)  - vector of lags for evaluation (seconds)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd2ccf.m 4095 2010-10-22 19:37:51Z karl $
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(csd)
    for i = 1:length(csd)
       [ccfi,pst] = spm_csd2ccf(csd{i},Hz);
       ccf{i}     = ccfi;
    end
    return
end
 
% Nyquist
%--------------------------------------------------------------------------
Ny    = 256;
g     = zeros(Ny,1);
Hz    = fix(Hz);
 
% Fourier transform cross-spectral density
%==========================================================================
for i = 1:size(csd,2)
    for j = 1:size(csd,3)
        g(Hz)      = csd(:,i,j);
        f          = ifft([0; g; flipud(conj(g))]);
        ccf(:,i,j) = real(fftshift(f));
    end
end
 
 
% Compute time bins if necessary
%--------------------------------------------------------------------------
N   = Ny;
pst = (-N:N)/N/2;
