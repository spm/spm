function [Y,FS] = spm_DEM_play_song(qU,T);
% displays the song-bird images specified by the states in qU
% FORMAT [Y,FS] = spm_DEM_play_song(qU,T);
%
% qU   - conditional moments of states (see spm_DEM)
% T    - number of seconds over which to play the sound
%
% Y    - sound image
% FS   - sampling rate (Hz)
%
% A button press on the spectrogram will play the song
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_play_song.m 6044 2014-06-14 10:22:46Z karl $
 
% load frequency modes
%--------------------------------------------------------------------------
try, T;           catch, T = 2;  end
try, v = qU.v{1}; catch, v = qU; end
 
% create sound image
%==========================================================================
[Nm m] = size(v);
 
% frequencies
%--------------------------------------------------------------------------
Hf  = 5000;                                % upper frequency (Hz)
Lf  = 2500;                                % lower frequency (Hz)
Nf  = 64;                                  % number of frequency bin
Hz  = linspace(Lf,Hf,64)';                 % frequencies
FS  = 2*Hz(end);                           % sampling rate (Hz)
k   = Hz/Hz(1);                            % cycles per window
n   = FS/Hz(1);                            % window length
N   = FS*T;                                % number of sonogram bins
R   = fix(N/m);                            % interpolation factor
N   = R*m;
pst = (1:N)/FS;                            % peristimulus time
sf  = 2*64^2;                              % dispersion of frequencies
 
 
% resample temporal modes
%--------------------------------------------------------------------------
for i = 1:Nm
    V(i,:) = spm_interp(v(i,:),R);
end
            
% create sonogram sound
%--------------------------------------------------------------------------
b     = V(1,:);                            % amplitude modulation
f     = V(2,:);                            % frequency modulation
b     = abs(b);
b     = b/max(b);
b     = tanh((b - 1/2)*6) + 1;
f     = 32*f + Lf;
S     = sparse(Nf,N);
for i = 1:N
    s      = b(i)*exp(-(Hz - f(i)).^2/sf);
    s      = sparse(s.*(s > exp(-4)));
    S(:,i) = s;
end


% inverse Fourier transform
%--------------------------------------------------------------------------
Y   = spm_iwft(S,k,n);
Y   = Y/max(Y);

 
% Graphics
%==========================================================================
colormap('pink')
imagesc(pst,Hz,abs(S))
axis xy
xlabel('time (sec)')
ylabel('Frequency (Hz)')
 
% set sound data
%--------------------------------------------------------------------------
h      = get(gca,'Children');
set(h(1),'Userdata',{Y,FS})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
