function [s] = spm_wavshow(I,WAV)
% rendering of a WAV sequence
% FORMAT [s] = spm_wavshow(I,WAV)
% I   - sonogram
% WAV - structure
% s   - sound
%
% This auxiliary routine plays the sound in a time frequency matrix (I),
% following an inverse wavelet transform. If called with no output
% arguments, it will also plot the sonogram (time frequency representation)
% and a sound file (timeseries) sampled from the sonogram (and sets the
% ButtonDownFcn of the sonogram image to play the sound).
%__________________________________________________________________________

% inverse wavelet transform
%--------------------------------------------------------------------------
s = spm_iwft(I,WAV.k,WAV.n);


if nargout, return, end

% image if no outputs 
%--------------------------------------------------------------------------
subplot(4,1,2)
imagesc(log(abs(I) + exp(-8))), xlabel('time (bins)'), ylabel('frequency')
title('Continuous Wavelet Transform','FontSize',12)

% Place movie in graphic subject
%--------------------------------------------------------------------------
h = get(gca,'Children');
set(h(1),'Userdata',{s,WAV.Fs})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% plot sound over seconds
%--------------------------------------------------------------------------
t = (1:numel(s))/WAV.Fs;
subplot(4,1,1)
plot(t,s), xlabel('time (seconds)'), ylabel('amplitude')
title('Sound','FontSize',12)
set(gca,'XLim',[t(1) t(end)])

% play sound
%--------------------------------------------------------------------------
soundsc(s,WAV.Fs);

return