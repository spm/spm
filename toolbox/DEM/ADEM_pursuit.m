% Slow pursuit under active inference: 
%__________________________________________________________________________
% This demo illustrates slow pursuit eye movements under active inference. 
% Its focus is on frames of references and the entrainment of gaze-
% direction by the motion of a visual target. The generative process (and 
% model) is based upon the itinerant trajectory of a target (in Cartesian 
% coordinates) produced by Lotka-Volterra dynamics. The agent expects its 
% sampling (in polar coordinates) to be centred on the target. Here, the 
% agent is equipped with a model of the trajectory and the oculomotor 
% plant. This means it represents both the location of the target and the 
% mapping from target location (in relation to a fixation point) to 
% egocentric polar coordinates. We simulate behavioural (saccadic) and
% electrophysiological (ERP) responses to expected and unexpected changes
% in the direction of a target moving on the unit circle.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_pursuit.m 4339 2011-06-02 13:04:49Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.a(:) - attractor (SHC) states
%
% v    - causal states
%   v(1) - not used
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%--------------------------------------------------------------------------
 
 
% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
n   = 8;                                      % number of attractors
a   = [1:n]*2*pi/n;                           % angles on unit circle
P   = [cos(a); sin(a)];
n   = size(P,2);                              % number of attractors
x.o = [0;0];                                  % oculomotor angle
x.x = [0;0];                                  % target location
x.a = sparse(1,1,4,n,1) - 6;                  % attractor (SHC) states
 
 
% Recognition model
%==========================================================================
M(1).E.s = 1;                                 % smoothness
M(1).E.n = 3;                                 % order of 
M(1).E.d = 2;                                 % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = 'spm_fx_dem_pursuit';              % plant dynamics
M(1).g   = 'spm_gx_dem_pursuit';              % prediction
 
M(1).x   = x;                                 % hidden states
M(1).V   = exp(6);                            % error precision
M(1).W   = exp(8);                            % error precision
M(1).pE  = P;
 
 
% level 2:
%--------------------------------------------------------------------------
M(2).v  = 0;                                  % inputs
M(2).V  = exp(2);
 
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_pursuit';
G(1).g  = 'spm_gx_adem_pursuit';
G(1).x  = x;                                   % hidden states
G(1).V  = exp(16);                             % error precision
G(1).W  = exp(16);                             % error precision
G(1).pE = P;
 
% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                  % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(8);
 
 
% generate and invert
%==========================================================================
N       = 128;                                % length of data sequence
D       = 0;
DEM.G   = G;
DEM.M   = M;
DEM.C   = [ones(1,N/2 + D) -ones(1,N/2 - D)];
DEM.U   = [ones(1,N/2) -ones(1,N/2)];
DEM     = spm_ADEM(DEM);
 
 
% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
 
% now repeat but delaying the reversal (unexpected)
%--------------------------------------------------------------------------
DUM     = DEM;
D       = 16;
DUM.C   = [ones(1,N/2 + D) -ones(1,N/2 - D)];
DUM     = spm_ADEM(DUM);

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_dem_pursuit_movie(DEM,0)
spm_dem_pursuit_movie(DUM,2)
subplot(2,2,3), title('Unexpected','FontSize',16)
subplot(2,2,4), title('Unexpected','FontSize',16)

  
% show saccadic and (simulated) ERP responses, time-locked to reversal
%==========================================================================
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%---------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
iE    = [D:(N - D)];
iU    = iE + D;
t     = [1:length(iE)]*8;
 
% true displacement from target (in intrinsic coordinates)
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(t,sum(DEM.pU.v{1}([3 4],iE).^2),'g'), hold on
plot(t,sum(DUM.pU.v{1}([3 4],iU).^2),'r'), hold off
xlabel('time (ms)','FontSize',12)
ylabel('angular (squared) distance from target','FontSize',12)
title('Saccadic (behavioural) responses','FontSize',16)
legend({'Expected'})
box off
 
subplot(2,1,2)
plot(t,DEM.qU.z{1}(:,iE),'g'), hold on
plot(t,DUM.qU.z{1}(:,iU),'r'), hold off
xlabel('time (ms)','FontSize',12)
ylabel('(sensory) prediction error','FontSize',12)
title('Electrophysiological responses','FontSize',16)
box off
