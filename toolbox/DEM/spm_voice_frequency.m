function [I] = spm_voice_frequency(Y,FS,F0)
% Segmentation of timeseries at fundamental frequency
% FORMAT [I] = spm_voice_frequency(Y,FS,F0)
%
% Y    - timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
%
% I    - intervals (time bins): mean(I) = DI = FS/F0
%
% This routine  identifies the the sampling intervals at the fundamental
% frequency, based upon the maxima after band-pass filtering around F0;
% namely, inflection or fluctuations in fundamental wavelength (i.e.,
% glottal pulse rate).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_frequency.m 7750 2019-12-05 17:54:29Z spm $


%% find fundamental frequencies
%==========================================================================
global VOX

% Fourier transform
%--------------------------------------------------------------------------
fY    = fft(Y(:));
nf    = length(fY);
w     = (1:nf)/(nf/FS);

% find fundamental frequency (f0 and fundamental epochs)
%--------------------------------------------------------------------------
R0    = F0/8;                                % standard deviation of F0(Hz)
bY    = fY.*exp(-(w(:) - F0).^2/(2*(R0)^2)); % bandpass filter (at F0)
[m,i] = max(abs(bY));                        % peak frequency
bY(i) = bY(i) + m/32;                        % add fundamental frequency
sY    = imag(ifft(bY));                      % filter timeseries

% intervals between maxima
%--------------------------------------------------------------------------
i     = find(diff(sY(1:end - 1)) > 0 & diff(sY(2:end)) < 0);
I     = diff(i)';


if ~isfield(VOX,'interval'), return, end

%% graphics
%==========================================================================
if VOX.interval
    
    % figure
    %----------------------------------------------------------------------
    spm_figure('GetWin','Voice (interval)'); clf;
    
    j    = 1:min(numel(sY),FS/2);
    b    = spm_zeros(sY);
    Y    =  Y/max(abs(Y(j)));
    sY   = sY/max(abs(sY(j)));
    b(i) = max(sY);
    pst  = 1000*j/FS;
    f0   = FS/mean(I);

    subplot(4,1,1)
    plot(pst,Y(j),pst,b(j),':r')
    str  = sprintf('Fundamental intervals (F0 = %.0f Hz)',f0);
    title(str,'FontSize',16), xlabel('time (ms)'), ylabel('amplitude')
    
    subplot(4,1,2)
    plot(pst,sY(j),pst,b(j),':r')
    str  = sprintf('Band pass filtering around %.0f Hz',F0);
    title(str,'FontSize',16), xlabel('time (ms)'), ylabel('amplitude')
    
    % graphics
    %--------------------------------------------------------------------------
    subplot(4,2,5)
    plot(i(2:end)/FS,FS./I),          hold on
    plot([0 i(end)/FS],[F0 F0],'r:'), hold on
    plot([0 i(end)/FS],[f0 f0],'r '), hold off
    title('Instantaneous frequency','FontSize',16)
    xlabel('time (seconds)'), ylabel('Hz')
    
    % graphics - F0
    %--------------------------------------------------------------------------
    subplot(4,2,6)
    i     = find(w > 64 & w < 300);
    plot(w(i),abs(fY(i))),                  hold on
    plot([F0 F0],[0 max(abs(fY(i)))],'r:'),  hold on
    plot([f0 f0],[0 max(abs(fY(i)))],'r '), hold off
    title('Fundamental frequency','FontSize',16)
    xlabel('frequency (Hz)'), ylabel('spectral energy'), drawnow
    
end
