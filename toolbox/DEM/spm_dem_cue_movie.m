function spm_dem_cue_movie(DEM)
% creates a movie of cued pointing
% FORMAT spm_dem_cue_movie(DEM)
%
% DEM - DEM structure from reaching simulations
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dem_cue_movie.m 4187 2011-02-01 20:13:57Z karl $


% evaluate location (targets)
%--------------------------------------------------------------------------
L    = DEM.pP.P{1}.x;

% evaluate positions (motor plant)
%--------------------------------------------------------------------------
x    = tan(DEM.pU.x{1}(1:2,:));                  % finger location
c    = DEM.pU.x{1}(3:end,:);                     % target contrast
c    = c - min(c(:));
c    = c/max(c(:));
fin  = imread('pointfinger.jpg');
N    = size(x,2);
n    = size(L,2);

% movie
%--------------------------------------------------------------------------
s     = 2;
for i = 1:N

    cla
    axis image ij
    hold on

    % finger
    %----------------------------------------------------------------------
    imagesc(([-1 0] + .68)*s + x(1,i),([-1 0] + .96)*s + x(2,i),fin);
    hold on

    % targets
    %----------------------------------------------------------------------
    for j = 1:n
        plot(L(1,j),L(2,j),'.','MarkerSize',64,'color',[c(j,i) (1 - c(j,i)) 0])
    end
    axis([-1 1 -1 1]*2)
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
