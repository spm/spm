function [I] = spm_voice_frequency(Y,FS,F0)
% segmentation of timeseries at fundamental frequency
% FORMAT [I] = spm_voice_frequency(Y,FS,F0)
%
% Y    - timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
%
% I    - intervals (time bins): mean(I) = DI = FS/F0
%
% This routine  identifies the fundamental formant frequency that captures
% the greatest variability in power upto 6 kHz. It then estimates the
% fundamental frequency using a Hilbert transform and returns the sampling
% intervals at the fundamental frequency; i.e., inflection or fluctuations
% in fundamental wavelength (i.e., glottal pulse rate).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_frequency.m 7575 2019-04-21 16:47:39Z karl $



%% find fundamental frequencies
%==========================================================================

% Fourier transform
%--------------------------------------------------------------------------
fY    = fft(Y(:));
nf    = length(fY);
w     = (1:nf)/(nf/FS);

% find fundamental frequency (f0 and fundamental epochs)
%--------------------------------------------------------------------------
R0    = F0/4;                                % standard deviation of F0(Hz)
bY    = fY.*exp(-(w(:) - F0).^2/(2*(R0)^2)); % bandpass filter (at F0)
sY    = real(ifft(bY));                      % filter timeseries
iY    = imag(hilbert(sY));                   % Hilbert transform
I     = diff(find(diff(iY < 0) > 0))';       % intervals (bins)

return

%% auxiliary graphics
%==========================================================================
subplot(4,1,1)
i    = cumsum(I);
j    = 1:min(numel(sY),FS/2);
b    = spm_zeros(sY);
b(i) = max(sY);
pst  = 1000*j/FS;
plot(pst,sY(j),pst,b(j),':')
title(sprintf('Fundamental intervals (F0 = %.0f Hz)',FS/mean(I)),'FontSize',16)
xlabel('time (seconds)')

% graphics
%--------------------------------------------------------------------------
subplot(4,2,3)
plot(i/FS,FS./I), hold on
plot([0 i(end)/FS],[F0 F0],'r'), hold off
title(sprintf('Instantaneous frequency',F0),'FontSize',16)
xlabel('time (seconds)')

% graphics - F0
%--------------------------------------------------------------------------
subplot(4,2,4)
i     = find(w > 64 & w < 300);
f0    = FS/mean(I);
plot(w(i),abs(fY(i))),                  hold on
plot([F0 F0],[0 max(abs(fY(i)))],'r'),  hold on
plot([f0 f0],[0 max(abs(fY(i)))],'r:'), hold off
title(sprintf('Fundamental frequency',F0),'FontSize',16)
xlabel('time (seconds)')










