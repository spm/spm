function [G,Hz] = spm_s2csd(s,Hz)
% Converts eigenspectrum to cross spectral density
% FORMAT [csd,Hz] = spm_s2csd(s,Hz)
%
% s                     - (DCM) eigenspectrum
% Hz   (n x 1)          - vector of frequencies (Hz)
%
% csd  (n,:,:)          - cross spectral density
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_s2csd.m 5013 2012-10-23 19:26:01Z karl $
 

% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*round(Hz(end)));
    N  = 1/dt;
    w  = round(linspace(Hz(1),Hz(end),length(Hz)));
catch
    N  = 128;
    dt = 1/N;
    w  = 1:N/2;
    Hz = w;
end
t      = (1:N)*dt;

% imaginary eigenmodes (s)
%--------------------------------------------------------------------------
s     = s(imag(s) > 0);
for i = 1:length(s)
    
    % transfer (S) functions (FFT of kernel)
    %----------------------------------------------------------------------
    S      = 1./(1j*(2*pi*w - imag(s(i))) - real(s(i))); 
    G(:,i) = S.*conj(S);
    
end

