function [I,DJ] = spm_voice_frequency(Y,FS,F0)
% segmentation of timeseries at fundamental frequency
% FORMAT [I,DJ] = spm_voice_frequency(Y,FS,F0)
%
% Y    - (spectral) timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
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
% $Id: spm_voice_frequency.m 7535 2019-03-03 20:27:09Z karl $



%% find fundamental frequencies
%==========================================================================

% Fourier transform
%--------------------------------------------------------------------------
fY    = fft(Y(:));
nf    = length(fY);
w     = (1:nf)/(nf/FS);

% find fundamental frequency (f0 and fundamental epochs)
%--------------------------------------------------------------------------
R0    = F0/8;                                % standard deviation of F0(Hz)
bY    = fY.*exp(-(w(:) - F0).^2/(2*(R0)^2)); % bandpass filter (at F0)
sY    = real(ifft(bY));                      % filter timeseries
iY    = imag(hilbert(sY));                   % Hilbert transform
I     = diff(find(diff(iY < 0) > 0))';       % intervals (bins)

if nargout < 2, return, end


%% find formant frequencies
%==========================================================================

% find fundamental formant frequency (F1) by maximising the standard
% deviation of power (absolute value) of discretely sampled frequencies
%--------------------------------------------------------------------------
N     = 64;                                  % number of harmonics
F1    = 128;                                 % formant frequency (Hz)
R1    = F1/16;                               % half range (Hz)
f1    = (-R1:R1) + F1;
sS    = spm_zeros(f1);
sF    = spm_conv(abs(fY),32*nf/FS);
for i = 1:length(f1)
    j     = round((1:N)*f1(i)*nf/FS);
    sS(i) = std(sF(j));                      % s.d. for this F1
end
[d,i] = max(sS);
F1    = f1(i);                               % formant frequency (Hz)
DJ    = FS/F1;                               % mean interval (bins)

return

%% graphics
%--------------------------------------------------------------------------
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


%% graphics
%--------------------------------------------------------------------------
subplot(4,1,3);
j    = round((1:N)*F1*nf/FS);
i    = 1:j(end);
b    = spm_zeros(sF);
b(j) = max(sF);
plot(w(i),abs(sF(i)),w(i),abs(b(i)),'r:')
title(sprintf('Formant frequencies (F1 = %.0f Hz)',F1),'FontSize',16)
xlabel('frequency (hertz)')

% graphics
%--------------------------------------------------------------------------
subplot(4,1,4);
plot(f1,sS)
title('Relative power','FontSize',16)
xlabel('formant frequency (hertz)')
ylabel('standard deviation')










