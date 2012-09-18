function [csd,Hz] = spm_s2csd(S,Hz)
% Converts eigenspectrum to cross spectral density
% FORMAT [csd,Hz] = spm_s2csd(S,Hz)
%
% S                     - (SPM) eigenspectrum
% Hz   (n x 1)          - vector of frequencies (Hz)
%
% csd  (n,:,:)          - cross spectral density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_s2csd.m 4936 2012-09-18 19:47:55Z karl $
 

% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*round(Hz(end)));
    N  = 1/dt;
    w  = round(linspace(Hz(1),M.Hz(end),length(Hz)));
catch
    N  = 128;
    dt = 1/N;
    w  = 1:N/2;
    Hz = w;
end
t  = (1:N)*dt;

% unpack specturm
%--------------------------------------------------------------------------
s     = [];
for i = 1:length(S)
    if imag(S(i))
        s(end + 1) = S(i);
        s(end + 1) = conj(S(i));
    else
        s(end + 1) = S(i);
    end
end

for i = 1:length(s)
    
    % transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    K     = exp(s(i)*t);

    % transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    csd    = fft(K);
    G(i,:) = csd(w + 1).*conj(csd(w + 1));
    
end

csd = real(sum(G)');

