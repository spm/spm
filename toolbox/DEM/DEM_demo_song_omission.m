% Demo for bird songs: In this example, we show that DEM can not only
% estimate the hidden states of an autonomous system but can also
% deconvolve dynamic changes in its control parameters.  We illustrate
% this using a slow Lorentz attractor to drive a fast one; both  showing
% deterministic chaos.  We endow the simulations with a little ethological
% validity by using the states of the fast Lorentz attractor as control
% variables in the syrinx of a song bird (usually these would control a van
% der Pol oscillator model). We will look at the true and inferred songs
% with and without the last part missing.  When sonograms are displayed the
% song can be played by a mouse click on the image. The final plots show
% simulated event related potential to show that there is a marked
% responses (prediction error) of the system when an expected ‘syllable’ is
% omitted. This demonstrates the implicit sequence-decoding of input
% streams, using uncontrollable state-space models
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_omission.m 4170 2011-01-24 18:37:42Z karl $
 
 
% Hierarchical non-linear generative model (dynamic & chaotic)
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
M(1).E.K = exp(-2);
 
% level 1
%--------------------------------------------------------------------------
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
M(1).f   = ' [-P(1) P(1) 0; (v(1) - 4 - x(3)) -1 0; x(2) 0 -P(2)]*x/16;';
M(1).fx  = '([-P(1) P(1) 0; (v(1) - 4 - x(3)) -1 0; x(2) 0 -P(2)] + [0 0 0; 0 0 -x(1); 0 x(1) 0])/16;';
M(1).fv  = '[0; x(2); 0]/16;';
M(1).g   = 'x([2 3])';
M(1).gx  = [0 1 0; 0 0 1];
M(1).x   = x;
M(1).pE  = P;
M(1).V   = exp(2);
M(1).W   = exp(8);
 
 
% level 2
%--------------------------------------------------------------------------
P        = [10; 8/3];
x        = [0.9; 0.8; 30];
M(2).f   = ' [-P(1) P(1) 0; (32 - x(3)) -1 0; x(2) 0 -P(2)]*x/128;';
M(2).fx  = '([-P(1) P(1) 0; (32 - x(3)) -1 0; x(2) 0 -P(2)] + [0 0 0; 0 0 -x(1); 0 x(1) 0])/128;';
M(2).fv  = sparse(3,0);
M(2).g   = 'x(3)';
M(2).gx  = [0 0 1];
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
 
DEM        = spm_DEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% now repeat with new stimulus
%==========================================================================
 
% suppress last syllable
%--------------------------------------------------------------------------
DEMa   = DEM;
DEMa.pU.v{1}(:,90:end) = 0;
DEMa.Y = DEMa.pU.v{1};
 
% deconvolve
%--------------------------------------------------------------------------
DEMa   = spm_DEM(DEMa);
spm_DEM_qU(DEMa.qU,DEMa.pU)
 
 
% show songs and prediction error (ERP)
%==========================================================================
spm_figure('Getwin','Graphics');
clf, colormap('pink')
 
% first stimulus
%--------------------------------------------------------------------------
subplot(3,2,1)
spm_DEM_play_song(DEM.pU,N*dt);
title('stimulus (sonogram)','Fontsize',18)
 
subplot(3,2,3)
spm_DEM_play_song(DEM.qU,N*dt);
title('percept','Fontsize',18)
 
subplot(3,2,5)
spm_DEM_EEG(DEM,dt,[1 2],1);
title('ERP (error)','Fontsize',18)
axis([1 N*dt*1000 -100 100])
 
% first stimulus
%--------------------------------------------------------------------------
subplot(3,2,2)
spm_DEM_play_song(DEMa.pU,N*dt);
title('without last syllable','Fontsize',18)
 
subplot(3,2,4)
spm_DEM_play_song(DEMa.qU,N*dt);
title('percept','Fontsize',18)
 
subplot(3,2,6)
spm_DEM_EEG(DEMa,dt,[1 2],1);
title('with omission','Fontsize',18)
axis([1 N*dt*1000 -100 100])
