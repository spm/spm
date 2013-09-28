function [ccf,pst] = spm_csd2ccf(csd,Hz,N)
% Converts cross spectral density to cross covariance function
% FORMAT [ccf,pst] = spm_csd2ccf(csd,Hz,N)
%
% csd  (n,:,:)          - cross spectral density (cf, mar.P)
% Hz   (n x 1)          - vector of frequencies (Hz)
% N                     - number of (positive) time bins (default: 256)
%
% ccf                   - cross covariance functions
% pst  (N,1)            - vector of lags for evaluation (seconds)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd2ccf.m 5660 2013-09-28 21:39:11Z karl $
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(csd)
    for i = 1:length(csd)
       [ccfi,pst] = spm_csd2ccf(csd{i},Hz);
       ccf{i}     = ccfi;
    end
    return
end
 
% unpack time bins (for time-frequency responses)
%--------------------------------------------------------------------------
if ndims(csd) == 4
    for i = 1:size(csd,1)
       [ccfi,pst]   = spm_csd2ccf(squeeze(csd(i,:,:,:)),Hz);
       ccf(i,:,:,:) = ccfi;
    end
    return
end


% Nyquist
%--------------------------------------------------------------------------
if nargin < 3
    N = 256;
end
N     = max(length(Hz),N);
g     = zeros(N,1);
dt    = 1/(Hz(2) - Hz(1));
Hz    = ceil(Hz*dt);
 
% Fourier transform cross-spectral density
%==========================================================================
for i = 1:size(csd,2)
    if ismatrix(csd)
        g(Hz)      = csd(:,i);
        f          = ifft([0; g; flipud(conj(g))]);
        ccf(:,i)   = real(fftshift(f))*N/dt;
    else
        for j = 1:size(csd,3)
            g(Hz)      = csd(:,i,j);
            f          = ifft([0; g; flipud(conj(g))]);
            ccf(:,i,j) = real(fftshift(f))*N/dt;
        end
    end
end
 
% Compute time bins
%--------------------------------------------------------------------------
pst = dt*(-N:N)/N/2;
