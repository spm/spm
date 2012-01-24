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
% $Id: ADEM_occulomotor_delays.m 4626 2012-01-24 20:55:59Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%--------------------------------------------------------------------------
 
 
% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
x.o = [0;0];                                  % oculomotor angle
x.x = [0;0];                                  % target location
 
 
% Recognition model
%==========================================================================
M(1).E.s = 1;                                 % smoothness
M(1).E.n = 4;                                 % order of 
M(1).E.d = 2;                                 % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = '[atan(x.x) - x.o; (v - x.x)/8]';  % extrinsic coordinates
M(1).g   = '[x.o; atan(x.x) - x.o]';          % intrinsic coordinate
M(1).x   = x;                                 % hidden states
M(1).V   = exp(8);                            % error precision
M(1).W   = exp(8);                            % error precision
 
% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                                  % inputs
M(2).V  = exp(0);
 
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = '[a - x.o/16; (v - x.x)/8]';        % extrinsic coordinates 
G(1).g  = '[x.o; atan(x.x) - x.o]';           % intrinsic coordinates
G(1).x  = x;                                  % hidden states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 64;                                  % length of data sequence
DEM.G   = G;
DEM.M   = M;
DEM.C   = [1; 1]*spm_phi(((1:N) - N/3)/2);
DEM.U   = zeros(2,N);
DEM     = spm_ADEM(DEM);
 
 
% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
 

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_dem_pursuit_movie(DEM,0)

