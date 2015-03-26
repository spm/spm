function [ccf,pst] = spm_csd2ccf(csd,Hz,dt)
% Converts cross spectral density to cross covariance function
% FORMAT [ccf,pst] = spm_csd2ccf(csd,Hz,dt)
%
% csd  (n,:,:)          - cross spectral density (cf, mar.P)
% Hz   (n x 1)          - vector of frequencies (Hz)
% dt                    - samping interval (default = 1/(2*Hz(end)))
%
% ccf                   - cross covariance functions
% pst  (N,1)            - vector of lags for evaluation (seconds)
%
% Note that because this scheme uses FFT one can only change dt.
%
% See also: 
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd2ccf.m 6395 2015-03-26 15:05:04Z adeel $

% Nyquist
%--------------------------------------------------------------------------
if nargin < 3, dt  = 1/(2*Hz(end)); end
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(csd)
    for i = 1:length(csd)
       [ccfi,pst] = spm_csd2ccf(csd{i},Hz,dt);
       ccf{i}     = ccfi;
    end
    return
end
 
% unpack time bins (for time-frequency responses)
%--------------------------------------------------------------------------
if ndims(csd) == 4
    for i = 1:size(csd,1)
       [ccfi,pst]   = spm_csd2ccf(squeeze(csd(i,:,:,:)),Hz,dt);
       ccf(i,:,:,:) = ccfi;
    end
    return
end


% indices for FFT
%--------------------------------------------------------------------------
dw    = Hz(2) - Hz(1);
Hz    = Hz/dw;
ns    = 1/dt;
N     = ceil(ns/2/dw);
gj    = find(Hz > 0 & Hz < (N + 1));
gi    = gj + ceil(Hz(1)) - 1;
g     = zeros(N,1);

% Fourier transform cross-spectral density
%==========================================================================
for i = 1:size(csd,2)
    if ismatrix(csd)
        g(gi)      = csd(gj,i);
        f          = ifft([0; g; flipud(conj(g))]);
        ccf(:,i)   = real(fftshift(f))*N*dw;
    else
        for j = 1:size(csd,3)
            g(gi)      = csd(gj,i,j);
            f          = ifft([0; g; flipud(conj(g))]);
            ccf(:,i,j) = real(fftshift(f))*N*dw;
        end
    end
end
 
% Compute time bins
%--------------------------------------------------------------------------
pst = dt*(-N:N);

