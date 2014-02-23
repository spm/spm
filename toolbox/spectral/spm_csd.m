function [csd,Hz] = spm_csd(Y,Hz,ns)
% Cross spectral density using Welch's method
% FORMAT [csd,Hz] = spm_csd(Y,Hz,ns)
%
% Y    (:,m)            - data
% Hz   (n x 1)          - vector of frequencies (Hz)
% ns                    - sampling frequency (default = 2*Hz(end))
%
% csd  (n,:,:)          - cross spectral density (cf, mar.P)
%
% See: cpsd.m and
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd.m 5892 2014-02-23 11:00:16Z karl $

% Nyquist
%--------------------------------------------------------------------------
if nargin < 3, ns  = 2*Hz(end); end
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(Y)
    for i = 1:length(csd)
       csd{i} = spm_csd(Y{i},Hz,dt);
    end
    return
end

% indices for FFT
%--------------------------------------------------------------------------
N     = round(ns/2);
f     = 0:N;
ci    = find(f >= Hz(1) & f <= Hz(end));

if ns > 8
    win = hanning(ns/4);  % EEG
else
    win = hanning(ns*64); % fMRI
end

% cross-spectral density
%==========================================================================
for i = 1:size(Y,2)
    for j = i:size(Y,2)
        c          = cpsd(Y(:,i),Y(:,j),win,[],2*N,ns);
        csd(:,i,j) = c(ci);
        csd(:,j,i) = conj(c(ci));
    end
end


