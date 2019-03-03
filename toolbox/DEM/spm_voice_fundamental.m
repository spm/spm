function [F0,F1] = spm_voice_fundamental(Y,FS)
% estimates fundamental and format frequencies
% FORMAT [F0,F1] = spm_voice_fundamental(Y,FS)
%
% Y    - timeseries
% FS   - sampling frequency
%
% F0   - fundamental frequency (glottal pulse rate)
% F1   - fundamental frequency (formant)
%
% This routine  identifies the fundamental and formant frequencies . The
% fundamental frequency is the lowest frequency (between 100 Hz and 300 Hz)
% that corresponds to the glottal pulse rate. The first format frequency is
% the average frequency (usually around 1000 Hz) that contains the spectral
% energy of the first formant.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_fundamental.m 7535 2019-03-03 20:27:09Z karl $



% find fundamental frequencies
%==========================================================================
i    = 1:min(FS*8,numel(Y));                 % analyse first 8 seconds

% Fourier transform
%--------------------------------------------------------------------------
fY    = abs(fft(Y(i)));                      % Fourier transform
nf    = length(fY);                          
w     = (1:nf)/(nf/FS);

% find fundamental frequency (F0)
%--------------------------------------------------------------------------
i     = find(w > 60 & w < 340);              % range of F0
sY    = spm_conv(fY(i),32);
j     = find(diff(diff(sY,1) > 0) < 0);
k     = find(sY(j) > max(sY(j))/2,1,'first');
F0    = w(i(1) + j(k));

if nargout < 2, return, end


% find fundamental formant frequency (F1)
%--------------------------------------------------------------------------
i     = find(w > 500 & w < 1500);            % range of F1
sF    = spm_conv(hamming(numel(i)).*fY(i),F0);
[d,j] = max(sF);
F1    = w(i(1) + j);

return

% graphics
%--------------------------------------------------------------------------
subplot(3,1,1);
i    = find(w < 512);
plot(w(i),fY(i),'r'),   hold on
plot([F0 F0],[0 max(fY(i))]), hold off
title('Fundamental frequency (F0))','FontSize',16)
xlabel('frequency (hertz)')

% graphics
%--------------------------------------------------------------------------
subplot(3,1,2);
i    = find(w < 8000);
plot(w(i),fY(i),'r'), hold on
for j = 1:8
    plot([F1 F1]*j,[0 max(fY(i))])
    for k = 1:8
        plot([F1 F1]*(j - 1 + (k - 1)/8),[0 max(fY(i))],':')
    end
end, hold off
title('Relative power','FontSize',16)
xlabel('formant frequency (hertz)')
ylabel('standard deviation')









