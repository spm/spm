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
% $Id: spm_voice_fundamental.m 7600 2019-06-01 09:30:30Z karl $


% find fundamental frequencies
%==========================================================================
global VOX, try VOX.formant; catch, VOX.formant = 0; end
i     = 1:min(FS*8,numel(Y));                   % analyse first 8 seconds

% Fourier transform
%--------------------------------------------------------------------------
Y     = Y(i) - spm_conv(Y(i),FS/64);            % high pass filter
fY    = abs(fft(Y));                            % Fourier transform
nf    = length(fY);
dw    = FS/nf;
nw    = 4;
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
I     = fix(FS/F0);
W     = (1:I)*FS/I;
for i = 1:(nf/I - 1)
    j      = (1:I) + (i - 1)*I;
    G(:,i) = abs(fft(Y(j)));
end
G     = var(G,1,2);
[~,j] = max(G);
F1    = W(j);

if ~VOX.formant, return, end

% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Frequencies'); clf

subplot(3,1,1)
i    = find(w < F0*(nw + 1));
j    = find(w > R0(1) & w < R0(end));
plot(w(i),fY(i),'r:',w(j),sY(j),'b'), hold on
plot([F0 F0],[0 max(fY(i))])
for j = 1:nw
    plot([F0 F0]*j,[0 max(fY(i))],'-.b')
end, hold off
title(sprintf('Fundamental frequency (F0 = %0.0f)',F0),'FontSize',16)
xlabel('frequency (hertz)'), box off

subplot(3,2,5)
plot(R0,S,'b'),                hold on
plot([F0 F0],[0 max(S)],'r'),  hold off
title('Relative energy (F0)','FontSize',16)
xlabel('F0 frequency (hertz)'), spm_axis tight, box off

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
title(sprintf('Formant frequency (F1 = %0.0f)',F1),'FontSize',16)
xlabel('formant frequency (Hz)')
ylabel('Energy'), box off

subplot(3,2,6)
j     = find(W  < 8000);
plot(W(j),G(j)),               hold on
plot([F1 F1],[0 max(G)],'r'),  hold off
title('Relative energy (F1)','FontSize',16)
xlabel('F1 frequency (hertz)'), spm_axis tight, box off






