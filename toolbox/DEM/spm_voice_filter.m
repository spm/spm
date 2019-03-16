function [G,F0] = spm_voice_filter(Y,FS,F1,F2)
% segmentation of timeseries at fundamental frequency
% FORMAT [G,F0] = spm_voice_filter(Y,FS)
%
% Y    - timeseries
% FS   - sampling frequency
% F1   - lower frequency bound [default: 512  Hz]
% F2   - upper frequency bound [default: 4096 Hz]
%
% G    - power at acoutic frequencies
% F0   - fundamental frequency
%
% This routine uses a simple wavelet decomposition (complex
% Gaussian wavelets) to assess the power frequency range (250 - 5000 Hz).
% This can be used to identify the onset of a word or fast modulations of
% spectral energy aat the fundamental frequency F0 (100 - 300Hz).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_filter.m 7545 2019-03-16 11:57:13Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 3; F1 = 512;  end
if nargin < 4; F2 = 8096; end


% find acoutic energy using spm_wft
%==========================================================================
F0    = 256;
k     = 1:2:round(F2/F0);                    % cycles per window                                % Hz
k     = k(k > F1/F0 & k < F2/F0);            % Acoustic range
n     = round(FS/(F0/2))*2;                  % window length (F0 Hz)
g     = abs(spm_wft(Y,k,n));                 % wavlet transform

% instantaneous power (250 - 5000 Hz)
%--------------------------------------------------------------------------
G     = sum(g)';

if nargout < 2, return, end
    
% find fundamental frequencies
%==========================================================================

% get maxium in the range of F0 (100 - 300Hz)
%--------------------------------------------------------------------------
fG    = abs(fft(G(:)));
nf    = length(fG);
w     = (1:nf)/(nf/FS);
i     = find(w > 64 & w < 300);
[d,j] = max(fG(i));
F0    = w(i(1) + j - 1);


return

% graphics
%==========================================================================
global VOX
if ~VOX.fundamental, return, end

% time-frequency analysis
%--------------------------------------------------------------------------  
subplot(3,1,1)
i    = 1:round(min(size(g,2),FS));
pst  = i*1000/FS;
imagesc(pst,k*F0,g(:,i))
title('time-frequency analysis','FontSize',16), xlabel('time (ms)')

% Energy
%--------------------------------------------------------------------------
subplot(3,1,2);
i    = 1:size(g,2);
pst  = i/FS;
plot(pst,G), title('Energy','FontSize',16), xlabel('time (s)')

% Fundamental frequency
%--------------------------------------------------------------------------
subplot(3,1,3)
i     = find(w > 64 & w < 300);
plot(w(i),fG(i)), hold on
plot([F0 F0],[0 max(fG(i))],'r'), hold off
title('Fundamental frequency','FontSize',16)
xlabel('time (seconds)')







