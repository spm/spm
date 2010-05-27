function [f]= spm_dem_reach_plot(DEM)
% plots the trajectory of a two-joint arm
% FORMAT [f]= spm_dem_reach_plot(DEM)
%
% DEM - DEM structure from reaching simulations
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_reach_plot.m 3901 2010-05-27 16:14:36Z karl $

% evaluate positions
%--------------------------------------------------------------------------
x    = DEM.pU.x{1};
v    = DEM.pU.v{2};
N    = size(x,2);

T    = [v(1,end); v(2,end)];                    % target location
J    = spm_dem_reach_x2J(x);                    % joint location
J{2} = J{1} + J{2};


% plot
%--------------------------------------------------------------------------
hold on
g     = [1 1 1]*.8;
for i = 1:N
    c = [1 .8 .7]*(1 - 0.2*i/N);
    plot([0 J{1}(1,i)],[0 J{1}(2,i)],'color',c)
    plot([J{1}(1,i) J{2}(1,i)],[J{1}(2,i) J{2}(2,i)],'color',c)
end
plot(J{2}(1,:),J{2}(2,:),'LineWidth',2,'color',c/2)
plot(J{2}(1,1),J{2}(2,1),'.r','MarkerSize',32)
plot(J{2}(1,end),J{2}(2,end),'.g','MarkerSize',32)
plot(T(1),T(2),'.r','MarkerSize',32)
hold off
axis image ij
axis([-0.5 1.5 0 2])