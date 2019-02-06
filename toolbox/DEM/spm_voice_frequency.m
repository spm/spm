function [I,DJ] = spm_voice_frequency(Y,FS,F0)
% segmentation of timeseries at fundamental frequency
% FORMAT [I,DJ] = spm_voice_frequency(Y,FS,F0)
%
% Y    - timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
% F1   - fundamental frequency (formant)
%
% I    - intervals (time bins): mean(I) = DI = FS/F0
% DJ   - FS/F1
%
% This routine  identifies the fundamental formant frequency that captures
% the greatest variability in power upto 6 kHz. It then estimates the
% fundamental frequency using a Hilbert transform and returns the sampling
% intervals at the fundamental frequency; i.e., inflection or fluctuations
% in fundamental wavelength (i.e., glottal pulse rate).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_frequency.m 7528 2019-02-06 19:19:49Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3, F0 = 100; end
F1    = 96;                                  % median formant frequency (Hz)
R1    = F1/6;                                % half range (Hz)
R0    = F1/6;                                % standard deviation of F0 (Hz)

% find fundamental frequencies
%==========================================================================

% find fundamental formant frequency (F1) by maximising the standard
% deviation of power (absolute value) of discretely sampled frequencies
%--------------------------------------------------------------------------
N     = 64;                                  % number of harmonics
fY    = fft(Y);                              % Fourier transform
nf    = length(fY);
sF    = spm_conv(abs(fY),32*nf/FS);
w     = (1:nf)/(nf/FS);
f1    = (-R1:R1) + F1;
sS    = spm_zeros(f1);
for i = 1:length(f1)
    j     = round((1:N)*f1(i)*nf/FS);
    sS(i) = std(sF(j));                      % s.d. for this F1
end
[d,i] = max(sS);
F1    = f1(i);                               % formant frequency (Hz)
DJ    = FS/F1;                               % mean interval (bins)

% find fundamental frequency (f0 and fundamental epochs)
%--------------------------------------------------------------------------
bY    = fY.*exp(-(w(:) - F0).^2/(2*(R0)^2)); % bandpass filter (at F0)
sY    = real(ifft(bY));                      % filter timeseries
iY    = imag(hilbert(sY));                   % Hilbert transform
I     = diff(find(diff(iY < 0) > 0))';       % intervals (bins)

return

% graphics
%--------------------------------------------------------------------------
subplot(3,1,1);
j    = round((1:N)*F1*nf/FS);
i    = 1:j(end);
b    = spm_zeros(sF);
b(j) = max(sF);
plot(w(i),abs(sF(i)),w(i),abs(b(i)),'r:')
title(sprintf('Formant frequencies (F1 = %.0f Hz)',F1),'FontSize',16)
xlabel('frequency (hertz)')

% graphics
%--------------------------------------------------------------------------
subplot(3,1,2);
plot(f1,sS)
title('Relative power','FontSize',16)
xlabel('formant frequency (hertz)')
ylabel('standard deviation')

% graphics
%--------------------------------------------------------------------------
subplot(3,1,3)
i    = cumsum(I);
j    = 1:3000;
b    = spm_zeros(Y);
b(i) = max(Y);
pst  = 1000*j/FS;
F0   = FS/mean(I);
plot(pst,Y(j),pst,b(j),':')
title(sprintf('Fundamental intervals (F0 = %.0f Hz)',F0),'FontSize',16)
xlabel('time (seconds)')









