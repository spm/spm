function spm_dem_search_trajectory(DEM)
% plots visual search in extrinsic and intrinsic coordinates
% FORMAT spm_dem_search_trajectory(DEM)
%
% DEM - {DEM} structures from visual search simulations
%
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception - x)
%   g(2) - oculomotor angle (proprioception - y)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_search_trajectory.m 4595 2011-12-19 13:06:22Z karl $
 
 
% Preliminaries
%--------------------------------------------------------------------------
clf, global STIM
N  = length(DEM);
S  = spm_read_vols(STIM.U);
 
% Stimulus
%======================================================================
Dx = STIM.U.dim(1)/2;
Dy = STIM.U.dim(2)/2;
a  = [];
q  = [];
c  = [];

subplot(2,2,1); hold off
image((S + 1)*32), axis image, hold on
 
for i = 1:N
    
    % i-th saccade - position
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.x{1}(1:2,:)*16;
    T  = length(pU);
    
    % eye movements in extrinsic coordinates
    %======================================================================
    subplot(2,2,1)
    plot(pU(2,:) + Dy,pU(1,:) + Dx,'r-','LineWidth',2)
    plot(pU(2,T) + Dy,pU(1,T) + Dx,'ro','MarkerSize',16)
    
    % Free energy
    %======================================================================
    F(i) = DEM{i}.F;
    
end
 
% Free energy
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(F)
title('Free-energy','FontSize',16)
xlabel('saccade')
axis square
