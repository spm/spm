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
% $Id: spm_csd2ccf.m 5837 2014-01-18 18:38:07Z karl $
 
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
if nargin < 3, dt  = 1/max(Hz)/2; end

ds    = Hz(2) - Hz(1);
N     = ceil(1/dt/2/ds);
gi    = ceil(Hz/ds);
gj    = 1:length(gi);
j     = gi > 1 & gi <= N;
gi    = gi(j);
gj    = gj(j);
g     = zeros(N,1);

% Fourier transform cross-spectral density
%==========================================================================
for i = 1:size(csd,2)
    if ismatrix(csd)
        g(gi)      = csd(gj,i);
        f          = ifft([0; g; flipud(conj(g))]);
        ccf(:,i)   = real(fftshift(f))*N*ds;
    else
        for j = 1:size(csd,3)
            g(gi)      = csd(gj,i,j);
            f          = ifft([0; g; flipud(conj(g))]);
            ccf(:,i,j) = real(fftshift(f))*N*ds;
        end
    end
end
 
% Compute time bins
%--------------------------------------------------------------------------
pst = dt*(-N:N);
