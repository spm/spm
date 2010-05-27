% This demo illustrates how action can fulfill prior expectations by 
% explaining away sensory prediction errors prescribed by desired movement 
% trajectories. In this example a two-joint arm follows a stable 
% heteroclinic channel, prescribed by a set of fixed point attractors. The 
% locations of the successive (autovitiated) attractors are defined by 
% parameters. The ensuing trajectories are illustrated here in terms of 
% synthetic writing.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_writing.m 3901 2010-05-27 16:14:36Z karl $


% hidden causes and states
%==========================================================================
% x    - hidden states
%   x.x(1) - joint angle
%   x.x(2) - joint angle
%   x.x(3) - angular velocity
%   x.x(4) - angular velocity
%
%   x.a(1) - attraction (location 1)
%   x.a(2) - attraction (location 2)
%   x.a(3) - attraction (location 3)
%    ...
%
% v    - causal states
%   v(1) - not used
%
%--------------------------------------------------------------------------
clear DEM

% parameters (locations) of trajectory
%--------------------------------------------------------------------------
P = [1.0  1.0;
     1.1  1.2;
     1.0  0.4;
     1.0  1.0;
     1.2  0.9;
     1.0  0.7;
     0.9  0.8;
     1.3  1.0]';
n = size(P,2);


% Recognition model (linear for expediency)
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of 
M(1).E.d = 2;                                 % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = 'spm_fx_dem_write';                % plant dynamics
M(1).g   = 'spm_gx_dem_write';                % prediction
 
M(1).x.x = [pi/2; pi/2; 0; 0];                % physical states
M(1).x.a = -[1:n]'/64;                       % attractor states
M(1).pE  = P;                                 % parameters of trajectory
M(1).V   = exp(8);                            % error precision
M(1).W   = exp(8);                            % error precision

 
% level 2: not used
%--------------------------------------------------------------------------
M(2).v  = [0];                                % inputs
M(2).V  = exp(8);
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_write';
G(1).g  = 'spm_gx_adem_write';
G(1).x  = [pi/2; pi/2; 0; 0;];                % physical states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 64;                                 % length of data sequence
DEM.G   = G;
DEM.M   = M;
DEM.C   = sparse(2,N);
DEM     = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% Graphics
%==========================================================================
spm_figure('GetWin','Graphics');
clf
 
subplot(2,1,1)
spm_dem_reach_plot(DEM)
title('trajectory','FontSize',16)
 
subplot(2,1,2)
spm_dem_reach_movie(DEM)
title('click on finger for movie','FontSize',16)
 

