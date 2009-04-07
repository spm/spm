% This demo illustrates how action can fulfil prior expectations by
% explaining away sensory prediction errors prescribed by desired movement
% trajectory. In this example a two-joint arm is trained to touch a traget
% so that sptronsous reaching occurs after training.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_reaching.m 3054 2009-04-07 19:22:49Z karl $

clear

% Recognition model (linear for expediency)
%==========================================================================
M(1).E.s      = 1/2;                          % smoothness
M(1).E.n      = 6;                            % order of 
M(1).E.d      = 2;                            % generlsied motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = 'spm_fx_dem_reach';                 % plant dynamics
M(1).g  = 'spm_gx_dem_reach';                 % prediction

M(1).x  = [pi/2; pi/2; 0; 0;];
M(1).V  = exp(8);                             % error precision
M(1).W  = exp(8);                             % error precision
 
% level 2: with noninformative priors on movement
%--------------------------------------------------------------------------
M(2).v  = [0; 0; 0];                          % inputs
M(2).V  = exp(0);
 
% generative model
%==========================================================================
G       = M;

% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_reach';
G(1).g  = 'spm_gx_adem_reach';
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0; 0];                          % inputs
G(2).a  = [0; 0];                             % action
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 64;                                 % length of data sequence
C       = sparse(3,N);
C(1,:)  = C(1,:) + .6;
C(2,:)  = C(2,:) + .6;
C(3,:)  = exp(-([1:N] - 32).^2/(8.^2));       % this is the prior cause

M(2).v  = C(:,1);


DEM.G   = G;
DEM.M   = M;
DEM.C   = C;
DEM.U   = sparse(3,N);
DEM0    = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM0.qU,DEM0.pU)
 

% Graphics
%==========================================================================
spm_figure('GetWin','Graphics');
clf

subplot(2,1,1)
spm_dem_reach_plot(DEM0)
title('trajectory','FontSize',16)

subplot(2,1,2)
spm_dem_reach_movie(DEM0)
title('click on finger for movie','FontSize',16)

return



% further simulations for paper
%==========================================================================









 
 
 