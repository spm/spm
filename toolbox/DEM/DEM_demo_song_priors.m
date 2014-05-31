function DEM_demo_song_priors
% Demo for a bird songs: In this example, we simulate local field potential
% using the prediction error from the song-bird example below. We look at
% these responses under natural stimuli and after removing the second
% level of the hierarchy to show it is necessary for veridical perception.
% We then repeat but omitting dynamical priors by forsaking generalised 
% coordinates
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_priors.m 6030 2014-05-31 13:09:24Z karl $
 
 
% hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================

% timing
%--------------------------------------------------------------------------
N        = 128;                      % length of stimulus (bins)
dt       = 1/64;                     % time bin (seconds)
 
% correlations
%--------------------------------------------------------------------------
M(1).E.s = 1;
M(1).E.K = exp(-2);


% level 1
%--------------------------------------------------------------------------
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
f        = '[-P(1) P(1) 0; (v(1) - 4 - x(3)) -1 0; x(2) 0 -P(2)]*x/16;';
M(1).f   = f;
M(1).g   = inline('x([2 3])','x','v','P');
M(1).x   = x;
M(1).pE  = P;
M(1).V   = exp(1);
M(1).W   = exp(8);
 
 
% level 2
%--------------------------------------------------------------------------
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
f        = '[-P(1) P(1) 0; (32 - x(3)) -1 0; x(2) 0 -P(2)]*x/128;';
M(2).f   = f;
M(2).g   = inline('x(3)','x','v','P');
M(2).x   = x;
M(2).pE  = P;
M(2).V   = exp(8);
M(2).W   = exp(8);
 
 
% create data
%==========================================================================
 
% create innovations & add causes
%--------------------------------------------------------------------------
DEM      = spm_DEM_generate(M,N);
 
% DEM estimation and display
%==========================================================================
DEM.M(1).x = [1; 1; 8];
DEM.M(2).x = [1; 1; 8];
 
% canoncial DEM
%--------------------------------------------------------------------------
DEMc   = DEM;
 
 
% without second level
%==========================================================================
DEMa   = DEM;
DEMa.M = DEMa.M(1);
 
% without generlised coordinates
%==========================================================================
DEMb   = DEM;
DEMb.M(1).E.n = 1;
 
% deconvolve
%--------------------------------------------------------------------------
DEMa   = spm_DEM(DEMa);
DEMb   = spm_DEM(DEMb);
DEMc   = spm_DEM(DEMc);


% show songs and prediction error (ERP)
%==========================================================================
spm_DEM_qU(DEMc.qU,DEMc.pU)
colormap('pink')

% Sonograms
%--------------------------------------------------------------------------
subplot(3,2,5)
spm_DEM_play_song(DEMc.pU ,N*dt);
title('simulus','Fontsize',18)
axis square

subplot(3,2,6)
spm_DEM_play_song(DEMc.qU ,N*dt);
title('percept','Fontsize',18)
axis square
drawnow


% Sonograms
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
clf, colormap('pink')

subplot(3,2,1)
spm_DEM_play_song(DEMc.qU ,N*dt);
title('percept (right click on image)','Fontsize',16)
 
subplot(3,2,3)
spm_DEM_play_song(DEMa.qU,N*dt);
title('no structural priors','Fontsize',16)
 
subplot(3,2,5)
spm_DEM_play_song(DEMb.qU,N*dt);
title('no dynamical priors','Fontsize',16)

 
% LFPs
%--------------------------------------------------------------------------
subplot(3,2,2)
spm_DEM_EEG(DEMc,dt,[1 2],1);
title('LFP','Fontsize',16)
 
subplot(3,2,4)
spm_DEM_EEG(DEMa,dt,[1 2],1);
title('LFP','Fontsize',16)
 
subplot(3,2,6)
spm_DEM_EEG(DEMb,dt,[1 2],1);
title('LFP','Fontsize',16)
drawnow, disp(' '),disp('Click sonograms to play songs'),disp(' ')

