function [Y,FS] = spm_DEM_play(qU,S,T);
% displays the sound images specified by the states in qU
% FORMAT [Y,FS] = spm_DEM_play(qU,S,T);
%
% qU   - conditional moments of states (see spm_DEM)
% S    - .mat file or structure 
%        U.U   - containing frequency modes (U) 
%        S.Hz  - and corresponding frequencies (FS)
% T    - number of second over which to play the sound
%
% Y    - sound image
% FS   - sampling rate (Hz)
%
% A button press on the spectrogram will play the sound
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_play.m 1380 2008-04-11 18:55:18Z karl $
 
% load frequency modes
%--------------------------------------------------------------------------
try
    T;
catch
    T = 2;
end
try
    S.U;
catch
    try
        S = load(file);
    catch
        S = load('BirdSong');
    end
end
 
% create sound image
%==========================================================================
v   = qU.v{1};
Hz  = S.Hz;
U   = S.U;
Nm  = size(U,2);


% frequencies
%--------------------------------------------------------------------------
FS  = 2*Hz(end);
k   = Hz/Hz(1);
n   = FS/Hz(1);
N   = FS*T;
R   = fix(N/size(v,2));
pst = [1:N]/FS;
 
% resample temporal modes
%--------------------------------------------------------------------------
for i = 1:Nm
    V(i,:) = interp(v(i,:),R);
    V(i,:) = V(i,:) - min(V(i,:));
end
            
% reconstituted sound
%--------------------------------------------------------------------------
SG  = U*V;
Y   = spm_iwft(SG,k,n);
Y   = Y/max(abs(Y));
 
% Graphics
%==========================================================================
imagesc(pst,Hz,abs(SG))
axis xy
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
 
% set sound data
%--------------------------------------------------------------------------
h      = get(gca,'Children');
set(h(1),'Userdata',{Y,FS})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
