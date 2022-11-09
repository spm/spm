function [csd,Hz] = spm_csd(Y,Hz,ns,PSD)
% Cross spectral density using Welch's method
% FORMAT [csd,Hz] = spm_csd(Y,Hz,ns)
%
% Y    (:,m)            - data
% Hz   (n x 1)          - vector of frequencies (Hz)
% ns                    - sampling frequency (default = 2*Hz(end))
% psd                   - 1 for power spectral density [default = 0]
%
% csd  (n,:,:)          - cross spectral density (cf, mar.P)
%
% See: cpsd.m and
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Nyquist
%--------------------------------------------------------------------------
if nargin < 3, ns  = 2*Hz(end); end
if nargin < 4, PSD = 0; end
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(Y)
    for i = 1:length(Y)
       Y{i} = spm_csd(Y{i},Hz,dt);
    end
    return
end

% indices for FFT
%--------------------------------------------------------------------------
Hz    = Hz(Hz <= ns/2);                   % Remove inestimable frequencies
dH    = Hz(2) - Hz(1);                    % frequency interval
f     = 0:dH:ns/2;                        % frequencies to estimate
N     = numel(f);                         % number of frequency points
ci    = find(f >= Hz(1) & f <= Hz(end));  % frequencies to retain

% cross-spectral density
%==========================================================================
for i = 1:size(Y,2)
    if PSD
        c        = pwelch(full(Y(:,i)),[],[],2*(N - 1),ns);
        csd(:,i) = c(ci);
    else
        for j = i:size(Y,2)
            c          = cpsd(full(Y(:,i)),full(Y(:,j)),[],[],2*(N - 1),ns);
            csd(:,i,j) = c(ci);
            csd(:,j,i) = conj(c(ci));
        end
    end
end


