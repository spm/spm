function ADEM_occlusion
% Slow pursuit and occlusion under active inference:
%__________________________________________________________________________
% This demo illustrates slow pursuit in the context of visual occlusion. We
% start with a simulation of canonical slow pursuit of a visual target 
% with sine wave motion. Crucially, the generative model is equipped with 
% a simple empirical prior encoding the hidden motion of the target (using 
% a linear oscillator, whose frequency is determined by a hidden cause). 
% We then examine failures of tracking and anticipation during occlusion 
% and when the target re-emerges from behind the occluder. We look at a 
% simulation in which the precision of the oscillator dynamics modelling 
% long-term behaviour of the target is reduced (cf., neuromodulatory 
% deficits in cortical areas encoding biological motion). This has a 
% an effect of producing a failure of pursuit, resulting in a catch-up
% saccade on target reappearance. The suppression of prior precision can
% however have beneficial effects when motion is itself unpredicted
% (as shown with differential pursuit performance under a reversal of 
% the trajectory towards the end of motion). Finally, we look at how prior 
% beliefs are acquired during exposure to the target – in terms of 
% cumulative inference on the hidden causes encoding the frequency of 
% periodic motion. This can be regarded as a high order form of evidence 
% accumulation. Importantly, this (experience-dependent) inference is
% markedly impaired by the simulated lesion to precision above. In other 
% words, a single failure of inference in modelling the motion of hidden 
% states can have secondary consequences – such as a failure to even 
% register and remember regularities. All these simulations are based upon 
% active inference; with the prior belief that the centre of gaze is 
% attracted to the same point responsible for target motion.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_occlusion.m 4804 2012-07-26 13:14:18Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor velocity
%   x.x(1) - target angle - extrinsic coordinates
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor velocity
%   g(:) - visual input - intrinsic coordinates
%--------------------------------------------------------------------------
 
 
% Set-up
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of
M(1).E.d = 1;                                 % generalised motion
 
 
% angular frequency of target motion
%--------------------------------------------------------------------------
w  = 2*pi/32;
 
 
% sensory mappings with and without occlusion
%--------------------------------------------------------------------------
g  = '[x.o; exp(-([-8:8]'' - x.x + x.o(1)).^2)*(x.x < 1/2)]';
h  = '[x.o; exp(-([-8:8]'' - x.x + x.o(1)).^2)]';
 
 
% oculomotor latencies (sinusoidal movement)
%==========================================================================
% Endow the model with internal dynamics (a simple oscillator) so that is
% recognises and remembers the trajectory to anticipate jumps in rectified
% sinusoidal motion. First, demonstrate canonical pursuit under occlusion:
 
% slow pursuit following with (second order) generative model
%--------------------------------------------------------------------------
x.o = [0;0];                                  % motor angle & velocity
x.x = 0;                                      % target location
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f = '[x.o(2); (v - x.o(1))/4 - x.o(2)/2; v - x.x]';
M(1).g = g;
M(1).x = x;                                   % hidden states
M(1).V = exp(4);                              % error precision
M(1).W = exp(4);                              % error precision
 
 
% level 2: With hidden (memory) states
%--------------------------------------------------------------------------
M(2).f  = '[x(2); -x(1)]*v/8';
M(2).g  = 'x(1)'; 
M(2).x  = [0; 0];                             % hidden states
M(2).V  = exp(4);                             % error precision
M(2).W  = exp(4);                             % error precision
 
% level 3: Encoding frequency of memory states (U)
%--------------------------------------------------------------------------
M(3).v = 0;
M(3).V = exp(4);
 
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f = '[x.o(2); a/4 - x.o(2)/8; v - x.x]';
G(1).g = g;
G(1).x = x;                                  % hidden states
G(1).V = exp(16);                            % error precision (errors)
G(1).W = exp(16);                            % error precision (motion)
G(1).U = sparse(1,[1 2],[1 1],1,19)*exp(4);  % motor gain
 
% second level
%--------------------------------------------------------------------------
G(2).v = 0;                                  % exogenous force
G(2).a = 0;                                  % action force
G(2).V = exp(16);
 
% Sine wave cause
%--------------------------------------------------------------------------
N      = 64;                                 % length of data sequence
dt     = 16;                                 % time step (ms)
t      = (1:N)*dt;                           % PST
 
DEM.M  = M;
DEM.G  = G;
DEM.C  = sin((1:N)*w).*((1:N) > 16);         % sinusoidal target motion
DEM.U  = zeros(1,N) + w*8;                   % prior beliefs
DEM    = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(3,2,1), title({'Slow pursuit:', 'prediction and error'},'FontSize',16)
subplot(3,2,2), title({'Occluded motion:', 'hidden states'},'FontSize',16)
 
% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_dem_occlusion_movie(DEM)
 
 
% repeat with simulated lesion to precision
%==========================================================================
SEM        = DEM;
SEM.M(2).W = exp(0);
SEM        = spm_ADEM(SEM);
spm_DEM_qU(SEM.qU,SEM.pU)
 
 
spm_figure('GetWin','Figure 3'); clf
spm_dem_occlusion_movie(SEM)
subplot(2,2,3), title({'Angular position:','reduced precision'},'FontSize',16)
subplot(2,2,4), title({'Angular velocity:','and catch-up saccade'},'FontSize',16)
 
 
% show improved tracking of unexpected trajectories under reduced precision
%==========================================================================
 
% remove occlusion and switch target trajectory after one cycle
%--------------------------------------------------------------------------
i          = 50;
DEM.M(1).g = h;
DEM.G(1).g = h;
DEM.C(i:N) = -DEM.C(i:N);
 
% reduce precisions and integrate
%--------------------------------------------------------------------------
SEM        = DEM;
SEM.M(2).W = exp(0);
 
DEM        = spm_ADEM(DEM);
SEM        = spm_ADEM(SEM);
 
spm_figure('GetWin','Figure 4'); clf
spm_dem_occlusion_movie(DEM)
subplot(2,2,3), hold on, subplot(2,2,4), hold on
spm_dem_occlusion_movie(SEM)
subplot(2,2,3), hold off, subplot(2,2,4), hold off
 
 
% illustrate inference on hidden cause (motion of target)
%==========================================================================
 
% allow for uncertainty about hidden causes (frequency of motion)
%--------------------------------------------------------------------------
DEM.M(3).V = exp(-4);
 
% remove occlusion
%--------------------------------------------------------------------------
DEM.M(1).g = h;
DEM.G(1).g = h;
 
% create a longer stimulus and reduce prior expectation
%--------------------------------------------------------------------------
N      = 128;
DEM.C  = sin((1:N)*w).*((1:N) > 16);
DEM.U  = zeros(1,N) + w/8;
DEM    = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 5');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(3,2,5), hold on, plot([1 N],[w w]*8,'-.k','LineWidth',4), hold off
 
 
% repeat with simulated lesion to precision
%==========================================================================
SEM        = DEM;
SEM.M(2).W = exp(0);
SEM        = spm_ADEM(SEM);
spm_DEM_qU(SEM.qU,SEM.pU)
 
spm_figure('GetWin','Figure 6'); clf
spm_DEM_qU(SEM.qU,SEM.pU)
subplot(3,2,4), title({'hidden states','with reduced precision'},'FontSize',16)
subplot(3,2,5), hold on, plot([1 N],[w w]*8,'-.k','LineWidth',4), hold off
axis([1 N -1/2 2])
