function [S] = DEM_birdsong(file)
% Create basis set for sounds
% FORMAT [S] = DEM_birdsong(file)
%
% file  - .wav file
%
% S.U   - h x 3 basis functions (Hz)
% S.V   - 3 x n basis functions (seconds)
% S.Hz  - s x 1 frequencies (Hz)
%
% Bird Song demo: These simple loads a .wav file of a real bird-song; and
% approximates the ensuing spectrogram with in terms of three
% time-frequency modes.  These modes are saved in BirdSong.mat (U) for
% illustrating DEM_demo_sequences
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% load bird song
%==========================================================================
if ~nargin
    file = fullfile(fileparts(mfilename('fullpath')),'lebi3.wav');
end
is_octave = strcmp(spm_check_version,'octave');
if is_octave || spm_check_version('matlab','8.0') > 0
    [Y,FS] = audioread(file);
else
    [Y,FS] = wavread(file); %#ok
end

Y   = Y(1:2^15);
T   = length(Y)/FS;                   % duration (seconds)

Nm  = 3;                              % number of frequency modes
cpw = 4;                              % minimum cycles per window
pst = (1:length(Y))/FS;               % peristimulus time (seconds)
Hz  = (1:64)*4000/64;                 % frequencies
k   = cpw*Hz/Hz(1);
n   = cpw*FS/Hz(1);

% windowed Fourier transform
%--------------------------------------------------------------------------
C   = spm_wft(Y(:),k,n);
C   = spm_conv(C,8,1);
Y   = spm_iwft(C,k,n);

% SVD
%--------------------------------------------------------------------------
U   = spm_svd(C*C');
U   = U(:,1:Nm);
V   = U'*C;
c   = U*V;

% reconstituted sound
%--------------------------------------------------------------------------
y   = spm_iwft(c,k,n);
y   = y/max(abs(y));

% output arguments
%--------------------------------------------------------------------------
S.U  = U;
S.V  = V;
S.Hz = Hz;


% Graphics
%==========================================================================
subplot(2,2,1)
imagesc(pst,Hz,abs(C))
axis xy
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('(with three modes')

% set sound data
%--------------------------------------------------------------------------
h      = get(gca,'Children');
set(h(1),'Userdata',{Y,FS})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

subplot(2,2,2)
imagesc(pst,Hz,abs(c))
axis xy
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('spectrogram')

% set sound data
%--------------------------------------------------------------------------
h      = get(gca,'Children');
set(h(1),'Userdata',{y,FS})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% plot modes
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(Hz,U)
xlabel('Frequency (Hz)')
title('Frequency modes')

subplot(2,2,4)
plot(pst,V)
xlabel('Time (sec)')
title('temporal modes')
