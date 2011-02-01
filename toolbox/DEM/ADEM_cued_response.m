% Cued responses under active inference: 
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
% egocentric polar coordinates.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_cued_response.m 4187 2011-02-01 20:13:57Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o  - oculomotor angle
%   x.x  - target locations (visual) - extrinsic coordinates (Cartesian)
%   x.a  - target contrast (attractiveness)
%
% v    - hidden causes
%
% P    - parameters
%
% g    - sensations:
%   g.o  - oculomotor angle (proprioception)
%   g.p  - target locations (visual) - intrinsic coordinates (polar)
%   g.c  - target contrast
%--------------------------------------------------------------------------

clear

% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
P.x = [-1  1 -1  1; -1  1  1 -1];

n   = size(P.x,2);                            % number of attractors

% hidden states (M)
%--------------------------------------------------------------------------
x.o = [0;0];                                  % oculomotor angle
x.x = [-1  1 -1  1; -1  1  1 -1];            % location states
x.a = sparse(1,1,4,n,1) - 6;                  % attractor (SHC) states

% hidden states (G)
%--------------------------------------------------------------------------
z.o = [0;0];                                  % oculomotor angle
z.a = sparse(1,1,4,n,1) - 6;                  % attractor (SHC) states


% Recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 2;                                 % order of 
M(1).E.d = 2;                                 % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = 'spm_fx_dem_cue';                  % plant dynamics
M(1).g   = 'spm_gx_dem_cue';                  % prediction
 
M(1).x   = x;                                 % hidden states
M(1).V   = exp(4);                            % error precision
M(1).W   = exp(8);                            % error precision
 
% level 2: not used
%--------------------------------------------------------------------------
M(2).f   = inline('spm_lotka_volterra(x,0)/8', 'x', 'v', 'P');
M(2).g   = inline('spm_softmax(x)/2',          'x', 'v', 'P');

M(2).x   = [1; 0];                            % hidden states
M(2).v   = 0;                                 % hidden causes
M(2).W   = exp(8);                            % error precision

% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_cue';
G(1).g  = 'spm_gx_adem_cue';
G(1).x  = z;                                  % hidden states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
G(1).pE = P;
 
% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                  % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 128;                                % length of data sequence
t       = [1:N]*8;
DEM.G   = G;
DEM.M   = M;
DEM.C   = sparse(1,N) + 1/2;
DEM     = spm_ADEM(DEM);

 
% dynamics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU)

subplot(3,2,5)
spm_dem_cue_movie(DEM)


 
 

