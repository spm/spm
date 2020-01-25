function [csd,Hz] = spm_ccf2csd(ccf,Hz)
% Converts cross covariance function to cross spectral density
% FORMAT [csd,Hz] = spm_ccf2csd(ccf,Hz)
%
% ccf  (N,:,:)          - cross covariance functions
% Hz   (n x 1)          - vector of frequencies (Hz)
%
% csd  (n,:,:)          - cross spectral density (cf, mar.P)
%
% See also: spm_???2???.m
%     ??? = {'ccf','csd','gew','mar','coh','mtf','ker','ssm','dcm'}
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ccf2csd.m 7774 2020-01-25 18:07:03Z karl $
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(ccf)
    for i = 1:length(ccf)
       csd{i}    = spm_ccf2csd(ccf{i},Hz);
    end
    return
end
 
% unpack time bins (for time-frequency responses)
%--------------------------------------------------------------------------
if ndims(ccf) == 4
    for i = 1:size(ccf,1)
       csd(i,:,:,:) = spm_ccf2csd(squeeze(ccf(i,:,:,:)),Hz);
    end
    return
end

% Frequencies
%--------------------------------------------------------------------------
ds    = Hz(2) - Hz(1);
gi    = ceil(Hz/ds) + 1;
N     = Hz(end);

% Fourier transform cross-spectral density
%==========================================================================
for i = 1:size(ccf,2)
    if ismatrix(ccf)
        g          = fft(ccf(:,i));
        csd(:,i)   = g(gi)/N/ds;
    else
        for j = 1:size(ccf,3)
            g          = fft(ccf(:,i));
            csd(:,i,j) = g(gi)/N/ds;
        end
    end
end
 
