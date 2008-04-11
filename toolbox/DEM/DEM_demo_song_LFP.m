% Demo for a bird songs: In this example, we simulate local field potential
% using the prediction error from the song-bird example above. We look at
% these responses under natural stimuli and after removing the second
% level of the hierarchy to show it is necessary for veridical perception.
% We then repeat but omitting dynamical priors by forsaking generalised 
% coordinates
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_LFP.m 1380 2008-04-11 18:55:18Z karl $
 
 
% hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================
spm_figure('Getwin','Graphics');
clear M
 
% timing
%--------------------------------------------------------------------------
N        = 128;                      % length of stimulus (bins)
dt       = 1/64;                     % time bin (seconds)
 
% correlations
%--------------------------------------------------------------------------
M(1).E.s = 1;
M(1).E.n = 6;
M(1).E.d = 2;
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
spm_DEM_qU(DEM.pU)
 
% DEM estimation and display
%==========================================================================
DEM.M(1).x = [0; 0; 8];
DEM.M(2).x = [1; 1; 8];
 
DEM    = spm_DEM(DEM);
 
% Repeat without second level
%==========================================================================
DEMa   = DEM;
DEMa.M = DEMa.M(1);
 
% deconvolve
%--------------------------------------------------------------------------
DEMa   = spm_DEM(DEMa);
 
% Repeat without generlised coordinates
%==========================================================================
DEMb   = DEM;
DEMb.M(1).E.n = 1;
 
% deconvolve
%--------------------------------------------------------------------------
DEMb   = spm_DEM(DEMb);
 
% show songs and predicti0n error (ERP)
%==========================================================================
spm_figure('Getwin','Graphics');
clf, colormap('pink')
 
% Sonograms
%--------------------------------------------------------------------------
subplot(3,2,1)
spm_DEM_play_song(DEM.qU ,N*dt);
title('natural stimulus','Fontsize',18)
 
subplot(3,2,3)
spm_DEM_play_song(DEMa.qU,N*dt);
title('no structural priors','Fontsize',18)
 
subplot(3,2,5)
spm_DEM_play_song(DEMb.qU,N*dt);
title('no dynamical priors','Fontsize',18)
 
% LFPs
%--------------------------------------------------------------------------
subplot(3,2,2)
spm_DEM_EEG(DEM ,dt);
title('LFP','Fontsize',18)
 
subplot(3,2,4)
spm_DEM_EEG(DEMa,dt);
title('LFP','Fontsize',18)
 
subplot(3,2,6)
spm_DEM_EEG(DEMb,dt);
title('LFP','Fontsize',18)
