function [G,w] = spm_s2csd(s,Hz)
% Convert eigenspectrum to cross spectral density
% FORMAT [csd,Hz] = spm_s2csd(s,Hz)
%
% s    (m x 1}        - eigenspectrum
% Hz   (n x 1)        - vector of frequencies (Hz)
%
% csd  (n,m)          - spectral density (of modes)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging
 

% frequencies of interest
%--------------------------------------------------------------------------
try
    w  = round(linspace(Hz(1),Hz(end),length(Hz)));
catch
    w  = 1:128;
end

% imaginary eigenmodes (s)
%--------------------------------------------------------------------------
for i = 1:length(s)
    
    % transfer (S) functions (FFT of kernel)
    %----------------------------------------------------------------------
    S      = 1./(1j*2*pi*w - s(i)); 
    G(:,i) = S.*conj(S);
    
end
