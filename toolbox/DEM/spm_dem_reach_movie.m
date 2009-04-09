function spm_dem_reach_movie(DEM)
% creates a movie of the trajectory of a two-joint arm
% FORMAT spm_dem_reach_movie(DEM)
%
% DEM - DEM structure from reaching simulations
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem_reach_movie.m 3058 2009-04-09 18:17:53Z karl $

% evaluate positions
%--------------------------------------------------------------------------
x    = DEM.pU.x{1};
v    = DEM.pU.v{2};

T    = [v(1,end); v(2,end)];                    % target location
J    = spm_dem_reach_x2J(x);                    % joint location
J{2} = J{1} + J{2};                             % finger location


fin  = imread('finger.jpg');
N    = size(x,2);

% movie
%--------------------------------------------------------------------------
s     = 1/4;

for i = 1:N
    
    cla
    axis image ij
    hold on
    imagesc(([-2 0] + 1/3)*s + J{2}(1,i),([-1 0] + 1/2)*s + J{2}(2,i),fin);
    plot([0 J{1}(1,i)],[0 J{1}(2,i)],'color',[1 .8 .7]*.8)
    plot([J{1}(1,i) J{2}(1,i)],[J{1}(2,i) J{2}(2,i)],'color',[1 .8 .7]*.8)
    plot(T(1),T(2),'.r','MarkerSize',32)
    axis([-0.5 1.5 0 2])
    hold off
    drawnow

    % save
    %----------------------------------------------------------------------
    M(i) = getframe(gca);
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
h = findobj(gca,'type','image');
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
