function [f]= spm_dem_reach_plot(DEM)
% plots the trajectory of a two-joint arm
% FORMAT [f]= spm_dem_reach_plot(DEM)
%
% DEM - DEM structure from reaching simulations
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_reach_plot.m 3054 2009-04-07 19:22:49Z karl $

% evaluate positions
%--------------------------------------------------------------------------
x    = DEM.pU.x{1};
v    = DEM.pU.v{2};

T    = [v(1,end); v(2,end)];                    % target location
J    = spm_dem_reach_x2J(x);                    % joint location
J{2} = J{1} + J{2};


% plot
%--------------------------------------------------------------------------
hold on
plot(J{2}(1,:),J{2}(2,:),'k')
plot(J{1}(1,:),J{1}(2,:),'k')
plot([0 J{1}(1,end)],[0 J{1}(2,end)],'color',[1 .8 .7]*.8)
plot([J{1}(1,end) J{2}(1,end)],[J{1}(2,end) J{2}(2,end)],'color',[1 .8 .7]*.8)
plot(J{2}(1,1),J{2}(2,1),'.r','MarkerSize',32)
plot(J{2}(1,end),J{2}(2,end),'.g','MarkerSize',32)
plot(T(1),T(2),'.r','MarkerSize',32)
hold off
axis image ij
axis([-0.5 1 0 1.5])