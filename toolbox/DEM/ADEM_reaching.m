function ADEM_reaching
% This demo illustrates how action can fulfil prior expectations by
% explaining away sensory prediction errors prescribed by desired movement
% trajectories. In this example a two-joint arm is trained to touch a target
% so that spontaneous reaching occurs after training.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_reaching.m 4804 2012-07-26 13:14:18Z karl $

% hidden causes and states
%==========================================================================
% x    - hidden states
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
% v    - causal states
%   v(1) - target location (x)
%   v(2) - target location (y)
%   v(3) - force (cue strength)
%--------------------------------------------------------------------------


% Recognition model (linear for expediency)
%==========================================================================
M         = struct;
M(1).E.s  = 1/2;                              % smoothness
M(1).E.n  = 4;                                % order of 
M(1).E.d  = 2;                                % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = 'spm_fx_dem_reach';                 % plant dynamics
M(1).g  = 'spm_gx_dem_reach';                 % prediction
 
M(1).x  = [pi/2; pi/2; 0; 0;];
M(1).V  = exp(8);                             % error precision
M(1).W  = exp(8);                             % error precision
 
% level 2: with non-informative priors on movement
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
G(1).U  = exp(8);                             % gain for action
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0; 0];                          % inputs
G(2).a  = [0; 0];                             % action
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 128;                                 % length of data sequence
C       = sparse(3,N);
C(1,:)  = C(1,:) + .6;                        % desired x
C(2,:)  = C(2,:) + .6;                        % desired y
C(3,:)  = exp(-((1:N) - 32).^2/(8.^2));       % cue strength
 
M(2).v  = C(:,1);
 
 
DEM.G   = G;
DEM.M   = M;
DEM.C   = C;
DEM.U   = sparse(3,N);
DEM     = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% Graphics
%==========================================================================
spm_figure('GetWin','Figure 1');
clf
 
subplot(2,1,1)
spm_dem_reach_plot(DEM)
title('trajectory','FontSize',16)
 
subplot(2,1,2)
spm_dem_reach_movie(DEM)
title('click on finger for movie','FontSize',16)
 
return
 
 
 
% further simulations for paper
%==========================================================================
 
% enable precision optimization
%--------------------------------------------------------------------------
Q{1}       = diag([1 1 0 0 0 0 0]);
Q{2}       = diag([0 0 0 0 0 1 1]);
M(1).E.nE  = 4;
M(1).E.nM  = 8;
M(1).V     = diag(exp([0 0 8 8 8 0 0]));               % error precision
M(1).Q     = Q;                                        % error components
 
 
% un-noisy
%--------------------------------------------------------------------------
G(1).V     = diag(exp([8 8 16 16 16 8 8]));            % error precision
M(1).hE    = [8 8];
DEM(1,1)   = DEM(1,1);
DEM(1,1).G = G;
DEM(1,1).M = M;
DEM(1,1)   = spm_ADEM(DEM(1,1));
 
% noisy position
%--------------------------------------------------------------------------
G(1).V     = diag(exp([8 8 16 16 16 4 4]));            % error precision
M(1).hE    = [8 4];
DEM(1,2)   = DEM(1,1);
DEM(1,2).G = G;
DEM(1,2).M = M;
DEM(1,2)   = spm_ADEM(DEM(1,2));
 
% noisy proprioception
%--------------------------------------------------------------------------
G(1).V     = diag(exp([4 4 16 16 16 8 8]));            % error precision
M(1).hE    = [4 8];
DEM(2,1)   = DEM(1,1);
DEM(2,1).G = G;
DEM(2,1).M = M;
DEM(2,1)   = spm_ADEM(DEM(2,1));
 
% noisy proprioception and position
%--------------------------------------------------------------------------
G(1).V     = diag(exp([4 4 16 16 16 4 4]));            % error precision
M(1).hE    = [4 4];
DEM(2,2)   = DEM(1,1);
DEM(2,2).G = G;
DEM(2,2).M = M;
DEM(2,2)   = spm_ADEM(DEM(2,2));
 
% show noisy proprioception
%--------------------------------------------------------------------------
spm_DEM_qU(DEM(2,1).qU,DEM(2,1).pU)
 
 
% overlay true values
%==========================================================================
spm_figure('GetWin','Figure 2');
clf
 
for i = 1:2
    for j = 1:2
        subplot(2,2,(j - 1)*2 + i)
        spm_dem_reach_plot(DEM(i,j))
    end
end
 
clf
for i = 1:2
    for j = 1:2
        subplot(4,4,(j - 1)*2 + i)
        bar(DEM(i,j).qH.h{1})
        set(gca,'YLim',[0 10])
    end
end
