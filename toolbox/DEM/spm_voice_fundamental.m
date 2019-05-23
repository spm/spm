function [F0,F1] = spm_voice_fundamental(Y,FS)
% Estimates and plots fundamental and format frequencies
% FORMAT [F0,F1] = spm_voice_fundamental(Y,FS)
%
% Y    - timeseries
% FS   - sampling frequency
%
% F0   - fundamental frequency (glottal pulse rate)
% F1   - fundamental frequency (formant)
%
% This auxiliary routine identifies the fundamental and formant
% frequencies. The fundamental frequency is the lowest frequency (between
% 100 Hz and 300 Hz) that corresponds to the glottal pulse rate. The first
% format frequency is the average frequency (usually around 1000 Hz) that
% contains the spectral energy of the first formant.
%
% This routine is not used for voice recognition but can be useful for
% diagnostic purposes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_fundamental.m 7597 2019-05-23 18:42:38Z karl $


% find fundamental frequencies
%==========================================================================
global VOX, try VOX.formant; catch, VOX.formant = 0; end
i     = 1:min(FS*8,numel(Y));                   % analyse first 8 seconds

% Fourier transform
%--------------------------------------------------------------------------
fY    = abs(fft(Y(i)));                         % Fourier transform
nf    = length(fY);
dw    = FS/nf;
nw    = 6;
sY    = spm_conv(fY,8/dw);
w     = (1:nf)*dw;

% find fundamental frequency (F0)
%--------------------------------------------------------------------------
R0    = 96:300;                                 % range of F0
for i = 1:numel(R0)
    DF   = R0(i)/dw;
    j    = (1:nw)*fix(DF);
    S(i) = sum(sY(j));
end
[j,i] = max(S);
F0    = R0(i);

if nargout == 1 && ~VOX.formant, return, end

% find fundamental formant frequency (F1)
%--------------------------------------------------------------------------
i     = find(w > 500 & w < 1500);               % range of F1
sF    = spm_conv(hamming(numel(i)).*fY(i),F0);
[~,j] = max(sF);
F1    = w(i(1) + j);

if nargout > 0 && ~VOX.formant, return, end

% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Frequencies');

subplot(3,1,1)
i    = find(w < F0*(nw + 1));
j    = find(w > R0(1) & w < R0(end));
plot(w(i),fY(i),'r:',w(j),sY(j),'b'), hold on
plot([F0 F0],[0 max(fY(i))])
for j = 1:nw
    plot([F0 F0]*j,[0 max(fY(i))],'-.b')
end, hold off
title(sprintf('Fundamental frequency (F0 = %0.0f)',F0),'FontSize',16)
xlabel('frequency (hertz)')

subplot(3,2,5)
plot(R0,S,'b'),                hold on
plot([F0 F0],[0 max(S)],'r'),  hold off
title('Relative energy (F0)','FontSize',16)
xlabel('F0 frequency (hertz)'), spm_axis tight

% graphics
%--------------------------------------------------------------------------
subplot(3,1,2)
i    = find(w < 8000);
plot(w(i),fY(i),'r'), hold on
for j = 1:8
    plot([F1 F1]*j,[0 max(fY(i))])
    for k = 1:8
        plot([F1 F1]*(j - 1 + (k - 1)/8),[0 max(fY(i))],':')
    end
end, hold off
title('Formant frequencies','FontSize',16)
xlabel('formant frequency (hertz)')
ylabel('Energy')









