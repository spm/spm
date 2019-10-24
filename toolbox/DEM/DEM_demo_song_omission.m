function DEM_demo_song_omission
% Demo for bird songs: In this example, we show that DEM can not only
% estimate the hidden states of an autonomous system but can also
% deconvolve dynamic changes in its control parameters.  We illustrate
% this using a slow Lorentz attractor to drive a fast one; both showing
% deterministic chaos.  We endow the simulations with a little ethological
% validity by using the states of the fast Lorentz attractor as control
% variables in the syrinx of a song bird (usually these would control a van
% der Pol oscillator model). We will look at the true and inferred songs
% with and without the last chirps missing.  The sonograms displayed
% can be played by a mouse click on the image.  Subsequent plots show
% simulated event-related potential to show that there is a marked
% responses (prediction error) of the system when an expected 'syllable' is
% omitted. This demonstrates the implicit sequence-decoding of input
% streams, using generative models based upon attractors.
% Having simulated normal omission-related responses, we then reduce the
% precision at the second level (on both hidden causes and states) and
% repeat the simulation. The result is an attenuation of the omission-
% related response or mismatch negativity. If we try to compensate by
% reducing the sensory precision, then the autonomous dynamics predicting
% the sequence of chirps supervenes, producing false inference. This
% can be thought of as a - crude - model of hallucinosis.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_omission.m 7679 2019-10-24 15:54:07Z spm $
 
 
% Hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================
 
% timing
%--------------------------------------------------------------------------
rng('default')
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
M(1).f   = @(x,v,P) [-P(1) P(1) 0; (v(1) - 4 - x(3)) -1 0; x(2) 0 -P(2)]*x/16;
M(1).g   = @(x,v,P) x([2 3]);
M(1).x   = x;
M(1).pE  = P;
M(1).V   = exp(2);
M(1).W   = exp(8);
 
 
% level 2
%--------------------------------------------------------------------------
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
M(2).f   = @(x,v,P) [-P(1) P(1) 0; (32 - x(3)) -1 0; x(2) 0 -P(2)]*x/128;
M(2).g   = @(x,v,P) x(3);
M(2).x   = x;
M(2).pE  = P;
M(2).V   = exp(8);
M(2).W   = exp(8);
 
 
% create data
%==========================================================================
 
% create innovations & add causes
%--------------------------------------------------------------------------
DEM      = spm_DEM_generate(M,N);
 
 
% reset initial hidden states and invert
%==========================================================================
DEM.M(1).x = [1; 1; 32];
DEM.M(2).x = [1; 1; 8];
DEM        = spm_DEM(DEM);
 
% illustrate responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
 
colormap('pink')
subplot(3,2,5)
spm_DEM_play_song(DEM.qU,N*dt);
title('percept','Fontsize',16)
 
 
% now repeat with new stimulus train, omitting the last chirps
%==========================================================================
dem   = DEM;
dem.pU.v{1}(:,86:end) = 0;
dem.Y = dem.pU.v{1};
 
% invert
%--------------------------------------------------------------------------
dem   = spm_DEM(dem);
 
spm_figure('GetWin','Figure 2');
spm_DEM_qU(dem.qU,dem.pU)
 
subplot(3,2,5)
spm_DEM_play_song(dem.qU,N*dt);
title('percept','Fontsize',16)
 
% show songs and prediction error (ERP)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
subplot(3,2,1)
spm_DEM_play_song(DEM.pU,N*dt);
title('stimulus (sonogram)','Fontsize',16)
 
subplot(3,2,3)
spm_DEM_play_song(DEM.qU,N*dt);
title('percept','Fontsize',16)
 
subplot(3,2,5)
spm_DEM_EEG(DEM,dt,[1 2]);
title('ERP (error)','Fontsize',16)
axis([1 N*dt*1000 -100 100])
 
% song with emissions
%--------------------------------------------------------------------------
subplot(3,2,2)
spm_DEM_play_song(dem.pU,N*dt);
title('without last syllable','Fontsize',16)
 
subplot(3,2,4)
spm_DEM_play_song(dem.qU,N*dt);
title('percept','Fontsize',16)
 
subplot(3,2,6)
spm_DEM_EEG(dem,dt,[1 2]);
title('with omission','Fontsize',16)
axis([1 N*dt*1000 -100 100])
 
 
% now repeat but with reduced precision to suppress violation responses
%==========================================================================
 
% illustrate previous omission-related responses and ERP
%--------------------------------------------------------------------------
dem.M(1).V = exp(2);
dem.M(1).W = exp(8);
dem.M(2).V = exp(16);
dem.M(2).W = exp(16);
dem        = spm_DEM(dem);

spm_figure('GetWin','Figure 4');
subplot(3,2,1)
spm_DEM_play_song(dem.qU,N*dt);
title('percept','Fontsize',16)
 
subplot(3,2,2)
spm_DEM_EEG(dem,dt,1);
title('response to violation','Fontsize',16)
axis([1 N*dt*1000 -100 100])
 
 
% reduce precision at first level (motion of hidden states)
%--------------------------------------------------------------------------
dem.M(1).V = exp(2);
dem.M(1).W = exp(8);
dem.M(2).V = exp(2);
dem.M(2).W = exp(2);
dem        = spm_DEM(dem);
 
spm_figure('GetWin','Figure 4');
subplot(3,2,3)
spm_DEM_play_song(dem.qU,N*dt);
title('percept','Fontsize',16)
 
subplot(3,2,4)
spm_DEM_EEG(dem,dt,1);
title('attenuated mismatch negativity','Fontsize',16)
axis([1 N*dt*1000 -100 100])
 
% and finally reduce sensory precision to compensate
%==========================================================================
dem.M(1).V = exp(-2);
dem.M(1).W = exp(8);
dem.M(2).V = exp(2);
dem.M(2).W = exp(2);
dem        = spm_DEM(dem);
 
spm_figure('GetWin','Figure 4');
subplot(3,2,5)
spm_DEM_play_song(dem.qU,N*dt);
title('percept','Fontsize',16)
 
subplot(3,2,6)
spm_DEM_EEG(dem,dt,1);
title('hallucination','Fontsize',16)
axis([1 N*dt*1000 -100 100])
 